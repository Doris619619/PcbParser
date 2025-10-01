#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <memory>
#include <cctype>
#include <algorithm>
#include <functional>
#include <unordered_set>

#include <iomanip>
#include <climits>

//哨兵
// === Tag helpers: [[FPID:x]] / [[SEID:x]] / [[VIID:x]] ===
static inline std::string makeTag(const char* kind, int id) {
    // kind = "FPID" / "SEID" / "VIID"
    return std::string("[[") + kind + ":" + std::to_string(id) + "]]";
}
static inline bool isTagged(const std::string& s, const char* kind) {
    const std::string p = std::string("[[") + kind + ":";
    return s.rfind(p, 0) == 0 && s.size() >= p.size() + 2 && s.substr(s.size() - 2) == "]]";
}



struct Node {
    std::string name;
    std::vector<std::string> parameters;
    std::vector<std::shared_ptr<Node>> children;

    Node(const std::string& n = "") : name(n) {}

    void addParameter(const std::string& param) {
        parameters.push_back(param);
    }

    void addChild(std::shared_ptr<Node> child) {
        children.push_back(child);
    }
};



// === 板层信息 ===
struct LayerInfo {
    int         id = -1;          // 例如 0, 31, 44 ...
    std::string name;             // 例如 F.Cu, B.Cu, Edge.Cuts ...
    std::string kind;             // 例如 signal, user
    std::string description;      // 例如 "Top", "Bottom", "B.Courtyard"；可能为空
};


class KiCadParser {
private:
    std::shared_ptr<Node> root;
    size_t pos;
    std::string content;
    int segmentCounter;
    int viaCounter;
    std::vector<std::shared_ptr<Node>> allSegments;
    std::vector<std::shared_ptr<Node>> allVias;
    std::vector<std::shared_ptr<Node>> allFootprints;
    std::vector<std::shared_ptr<Node>> allNodes;  // 用来存储所有解析的节点
    std::vector<LayerInfo> boardLayers;  // 顶层 (layers ...) 的解析结果

    // ===== 删除辅助：按 tstamp 精确匹配 =====
// 仅删除 name==expectName 且其子节点包含 (tstamp <exact>) 的节点
    size_t removeNodesByTstampRec(const std::shared_ptr<Node>& cur,
        const std::string& tstamp,
        const char* expectName) {
        if (!cur) return 0;
        size_t removed = 0;
        auto& ch = cur->children;
        const size_t before = ch.size();
        ch.erase(std::remove_if(ch.begin(), ch.end(),
            [&](const std::shared_ptr<Node>& c) {
                if (!c) return false;
                if (expectName && c->name != expectName) return false;
                for (const auto& gc : c->children) {
                    if (gc && gc->name == "tstamp" && !gc->parameters.empty() && gc->parameters[0] == tstamp)
                        return true;
                }
                return false;
            }), ch.end());
        removed += (before - ch.size());
        for (auto& c : ch) removed += removeNodesByTstampRec(c, tstamp, expectName);
        return removed;
    }


    // 解析顶层 (layers ...)；只扫 root 的直接子节点，避免与 pad 内的 (layers ...) 混淆
    void parseBoardLayers() {
        boardLayers.clear();
        if (!root) return;

        // 只在 root->children 里找名字为 "layers" 且有子节点的块（顶层板层表就是这种形态）
        for (const auto& ch : root->children) {
            if (!ch) continue;
            if (ch->name != "layers") continue;
            if (ch->children.empty()) continue; // pad 的 (layers "F.Cu" ...) 没有子节点

            for (const auto& lay : ch->children) {
                if (!lay) continue;

                // 子项形如： (0 "F.Cu" signal "Top")
                // 解析规则：name=编号，parameters[0]=层名，parameters[1]=类型，parameters[2]=可选描述
                int id = -1;
                try { id = std::stoi(lay->name); }
                catch (...) { /* 非数字就跳过 */ continue; }

                LayerInfo info;
                info.id = id;
                if (!lay->parameters.empty())      info.name = lay->parameters[0];
                if (lay->parameters.size() >= 2)   info.kind = lay->parameters[1];
                if (lay->parameters.size() >= 3)   info.description = lay->parameters[2];

                boardLayers.push_back(std::move(info));
            }
            break; // 顶层 (layers ...) 只会有一个，找到就停
        }
    }



    //添加via部分
    // —— 数字转字符串：尽量还原“人能读”的样子（不科学计数法、合适精度）
    static std::string fmtNum(double v) {
        std::ostringstream oss;
        oss << std::setprecision(15) << std::noshowpoint << std::defaultfloat << v;
        return oss.str();
    }

    // —— 从 [[KIND:123]] 里提取 id，提取失败返回 -1
    static int extractTagId(const std::string& p, const char* kind) {
        if (!isTagged(p, kind)) return -1;
        const size_t colon = p.find(':');
        const size_t end = p.rfind("]]");
        if (colon == std::string::npos || end == std::string::npos || end <= colon + 1) return -1;
        try { return std::stoi(p.substr(colon + 1, end - colon - 1)); }
        catch (...) { return -1; }
    }

    // —— 计算下一枚 via 的编号：现有 [[VIID:x]] 的最大值 + 1；若没有，则从 1 开始
    int nextViaId() const {
        int mx = 0;
        for (const auto& v : allVias) {
            if (!v || v->parameters.empty()) continue;
            int id = extractTagId(v->parameters.front(), "VIID");
            if (id > mx) mx = id;
        }
        return mx + 1;
    }


    // === 删除相关：匹配与递归移除 ===

// 判断节点的首参数是否是期望的精确哨兵，比如 [[VIID:12]]
    static inline bool hasExactTag(const std::shared_ptr<Node>& n, const char* kind, int id) {
        if (!n || n->parameters.empty()) return false;
        return n->parameters.front() == makeTag(kind, id);
    }

    // 递归从整棵树中删除：所有满足 (expectName && 首参数 == [[KIND:id]]) 的节点
    // 返回删掉的数量（可能是 0 / 1 / 多个）
    size_t removeNodesByTagRec(const std::shared_ptr<Node>& cur,
        const char* kind,
        int id,
        const char* expectName) {
        if (!cur) return 0;
        size_t removed = 0;

        // 先在当前节点的直接 children 里删
        auto& ch = cur->children;
        const size_t before = ch.size();
        ch.erase(std::remove_if(ch.begin(), ch.end(),
            [&](const std::shared_ptr<Node>& c) {
                if (!c) return false;
                if (expectName && c->name != expectName) return false;
                return hasExactTag(c, kind, id);
            }),
            ch.end());
        removed += (before - ch.size());

        // 再递归处理剩余孩子
        for (auto& c : ch) {
            removed += removeNodesByTagRec(c, kind, id, expectName);
        }
        return removed;
    }

    // 删除后不重新编号，只“重建缓存”列表，保持其它编号不变
    void rebuildCaches() {
        allNodes.clear();
        collectAllNodes(root);

        allVias.clear();
        findNodesRecursive(root, "via", allVias);

        allSegments.clear();
        findNodesRecursive(root, "segment", allSegments);

        allFootprints.clear();
        findNodesRecursive(root, "footprint", allFootprints);
    }


    // ----- 序列化辅助：判断数字 -----
    static bool isNumberToken(const std::string& s) {
        if (s.empty()) return false;
        bool seenDigit = false, seenDot = false, seenExp = false;
        size_t i = 0; if (s[i] == '+' || s[i] == '-') ++i;
        for (; i < s.size(); ++i) {
            char c = s[i];
            if (std::isdigit((unsigned char)c)) { seenDigit = true; continue; }
            if (c == '.') { if (seenDot || seenExp) return false; seenDot = true; continue; }
            if (c == 'e' || c == 'E') {
                if (seenExp || !seenDigit) return false;
                seenExp = true; seenDigit = false;
                if (i + 1 < s.size() && (s[i + 1] == '+' || s[i + 1] == '-')) ++i;
                continue;
            }
            return false;
        }
        return seenDigit;
    }

    // ----- 转义到双引号 -----
    static std::string escapeForQuotes(const std::string& s) {
        std::string out; out.reserve(s.size() + 8);
        for (char c : s) { if (c == '"' || c == '\\') out.push_back('\\'); out.push_back(c); }
        return out;
    }

    // ----- 是否需要引号 -----
    // 规则：空串 → 必须引号；数字 → 不引号；
    // 只含 [a-z0-9_\-./+:] → 不引号（全部小写/安全字符）；
    // 含大写或空白/其他符号 → 引号。
    static bool shouldQuote(const std::string& s) {
        if (s.empty()) return true;           // 关键：输出 ""，用于 group "" / net_name "" / outputdirectory ""
        if (isNumberToken(s)) return false;
        bool hasUpper = false;
        for (unsigned char uc : s) {
            char c = (char)uc;
            if (std::isupper(uc)) hasUpper = true;
            if (!(std::islower(uc) || std::isdigit(uc) ||
                c == '_' || c == '-' || c == '.' || c == '/' || c == '+' || c == ':')) {
                return true; // 有空白、引号、反斜杠或其他符号
            }
        }
        return hasUpper; // 有大写也加引号（如 "B.Cu", "Top"）
    }

    // ----- 参数渲染 -----
    static std::string renderParam(const std::string& p) {
        if (p.empty()) return "\"\"";                      // ★ 空字符串：强制输出 ""
        if (!shouldQuote(p)) return p;
        return std::string("\"") + escapeForQuotes(p) + "\"";
    }

    // ----- 递归写节点（处理 hide 的位置） -----
    void writeNodeRec(std::ostream& os,
        const std::shared_ptr<Node>& n,
        int indent,
        int indentStep) {
        if (!n) return;

        std::string indentStr((size_t)indent, ' ');

        // 1) 先把参数里的 "hide" 拆出来
        std::vector<std::string> paramsNoHide;
        int hideCount = 0;
        for (const auto& p : n->parameters) {
            if (p == "hide") ++hideCount;
            else paramsNoHide.push_back(p);
        }

        // 2) 写“(name + 普通参数)”
        os << indentStr << "(" << n->name;
        for (const auto& p : paramsNoHide)
            os << " " << renderParam(p);

        // 3) 如果没有子节点：把 hide 仍放在行尾，保持兼容
        if (n->children.empty()) {
            while (hideCount-- > 0) os << " hide";
            os << ")\n";
            return;
        }

        // 4) 有子节点：换行，然后在写完第一个子节点后插入一行 "hide"
        os << "\n";

        bool injectedHide = false;
        for (size_t i = 0; i < n->children.size(); ++i) {
            const auto& ch = n->children[i];
            writeNodeRec(os, ch, indent + indentStep, indentStep);

            if (!injectedHide && hideCount > 0) {
                std::string ind2((size_t)(indent + indentStep), ' ');
                while (hideCount-- > 0) {
                    os << ind2 << "hide\n";  // 紧跟在第一个“括号子项”之后
                }
                injectedHide = true;
            }
        }

        os << indentStr << ")\n";
    }


    

    


    //哨兵
    // --- 在 KiCadParser::private: 里 ---
    void stripAllTempIdsRec(const std::shared_ptr<Node>& n) {
        if (!n) return;
        if (!n->parameters.empty()) {
            const std::string& p0 = n->parameters.front();
            if (isTagged(p0, "FPID") || isTagged(p0, "SEID") || isTagged(p0, "VIID")) {
                const_cast<std::vector<std::string>&>(n->parameters).erase(n->parameters.begin());
            }
        }
        for (auto& ch : n->children) stripAllTempIdsRec(ch);
    }

    


    /*
    void skipWhitespace() {
        while (pos < content.size() && std::isspace(content[pos])) {
            pos++;
        }
    }
    */

    // 1) 建议把 isspace 的参数都转成 unsigned char，防止中文/高位字节引发UB
    void skipWhitespace() {
        while (pos < content.size() && std::isspace(static_cast<unsigned char>(content[pos]))) {
            ++pos;
        }
    }

    // 2) 新增：从当前位置读取一个“裸 token”，不跳过开头空白，不吞括号
    std::string readBareToken() {
        std::string token;
        while (pos < content.size()) {
            char c = currentChar();
            if (std::isspace(static_cast<unsigned char>(c)) || c == '(' || c == ')') break;
            token += c;
            ++pos;
        }
        return token;
    }

    char currentChar() {
        if (pos < content.size()) {
            return content[pos];
        }
        return '\0';
    }

    char nextChar() {
        pos++;
        return currentChar();
    }

    std::string readQuotedString() {
        if (currentChar() != '"') return "";
        ++pos;                       // 吃掉开头的 "
        std::string out;
        while (pos < content.size() && currentChar() != '"') {
            if (currentChar() == '\\') ++pos;     // 处理转义
            if (pos < content.size()) {
                out.push_back(currentChar());
                ++pos;
            }
        }
        if (currentChar() == '"') ++pos;          // 吃掉收尾的 "
        return out;                               // 可能是空串 ""
    }



    std::string readToken() {
        skipWhitespace();

        if (pos >= content.size()) {
            return "";
        }

        char c = currentChar();
        if (c == '(' || c == ')') {
            pos++;
            return std::string(1, c);
        }

        if (c == '"') {
            return "\"" + readQuotedString() + "\"";
        }

        std::string token;
        while (pos < content.size()) {
            c = currentChar();
            if (std::isspace(c) || c == '(' || c == ')') {
                break;
            }
            token += c;
            pos++;
        }

        return token;
    }

    std::shared_ptr<Node> parseNode() {
        skipWhitespace();

        if (currentChar() != '(') {
            return nullptr;
        }

        ++pos;
        skipWhitespace();

        std::string name;
        while (pos < content.size()) {
            char c = currentChar();
            if (std::isspace(static_cast<unsigned char>(c)) || c == '(' || c == ')') break;
            name += c;
            ++pos;
        }
        if (name.empty()) return nullptr;

        auto node = std::make_shared<Node>(name);

        while (pos < content.size()) {
            skipWhitespace();

            if (currentChar() == ')') { ++pos; break; }

            if (currentChar() == '(') {
                auto child = parseNode();
                if (child) node->addChild(child);
            }
            else {
                if (currentChar() == '"') {
                    // ★ 保留空串参数
                    std::string q = readQuotedString();
                    node->addParameter(q);
                }
                else {
                    std::string tok = readBareToken();   // 若没有这个函数就用 readToken()
                    if (!tok.empty()) node->addParameter(tok);
                }
            }
        }

        return node;
    }


    /*
    void findNodesRecursive(const std::shared_ptr<Node>& node,
        const std::string& targetName,
        std::vector<std::shared_ptr<Node>>& results) {
        if (!node) return;

        if (node->name == targetName) {
            results.push_back(node);
        }

        for (const auto& child : node->children) {
            findNodesRecursive(child, targetName, results);
        }
    }
    */
    /*
    void findNodesRecursive(const std::shared_ptr<Node>& node,
        const std::string& targetName,
        std::vector<std::shared_ptr<Node>>& results) {
        if (!node) return;

        // 将所有节点添加到 allNodes 中
        allNodes.push_back(node);

        // 查找特定名称的节点并将其添加到结果列表中
        if (node->name == targetName) {
            results.push_back(node);
        }

        // 递归遍历子节点
        for (const auto& child : node->children) {
            findNodesRecursive(child, targetName, results);
        }
    }
    */

    // 只做“按名查找”，不再往 allNodes 里塞任何东西
    void findNodesRecursive(const std::shared_ptr<Node>& node,
        const std::string& targetName,
        std::vector<std::shared_ptr<Node>>& results) {
        if (!node) return;
        if (node->name == targetName) {
            results.push_back(node);
        }
        for (const auto& child : node->children) {
            findNodesRecursive(child, targetName, results);
        }
    }

    // 专门用来收集所有节点
    void collectAllNodes(const std::shared_ptr<Node>& node) {
        if (!node) return;
        allNodes.push_back(node);
        for (const auto& child : node->children) {
            collectAllNodes(child);
        }
    }




    //给所有segment编号
    void findAndNumberSegments() {
        segmentCounter = 0;
        allSegments.clear();
        findNodesRecursive(root, "segment", allSegments);

        for (auto& segment : allSegments) {
            segmentCounter++;
            std::string tag = makeTag("SEID", segmentCounter);
            if (segment->parameters.empty() || !isTagged(segment->parameters.front(), "SEID")) {
                segment->parameters.insert(segment->parameters.begin(), tag);
            }
            else {
                segment->parameters.front() = tag; // 已有哨兵则覆盖
            }

        }
    }



    void findAndNumberVias() {
        viaCounter = 0;
        allVias.clear();
        findNodesRecursive(root, "via", allVias);
        for (auto& via : allVias) {
            ++viaCounter;
            std::string tag = makeTag("VIID", viaCounter);
            if (via->parameters.empty() || !isTagged(via->parameters.front(), "VIID")) {
                via->parameters.insert(via->parameters.begin(), tag);
            }
            else {
                via->parameters.front() = tag;
            }

        }
    }


    void findAndNumberFootprints() {
        int footprintCounter = 0;  // 用来给footprint编号
        allFootprints.clear();  // 清空原有的footprints
        findNodesRecursive(root, "footprint", allFootprints);  // 查找所有footprint节点

        for (auto& footprint : allFootprints) {
            footprintCounter++;
            std::string tag = makeTag("FPID", footprintCounter);
            if (footprint->parameters.empty() || !isTagged(footprint->parameters.front(), "FPID")) {
                footprint->parameters.insert(footprint->parameters.begin(), tag);
            }
            else {
                footprint->parameters.front() = tag;
            }

        }
    }


public:
    KiCadParser() : pos(0), segmentCounter(0), viaCounter(0) {}

    // ===== 导出 .kicad_pcb =====
// 默认保存前清除 [[FPID:]] / [[SEID:]] / [[VIID:]] 哨兵，避免污染文件
    bool saveAsKicadPcb(const std::string& outPath,
        bool stripTempIds = true,
        int indentStep = 2) {
        if (!root) return false;
        if (stripTempIds) stripAllTempIds();
        std::ofstream ofs(outPath, std::ios::binary);
        if (!ofs.is_open()) {
            std::cerr << "无法写出文件: " << outPath << "\n";
            return false;
        }
        writeNodeRec(ofs, root, /*indent=*/0, indentStep);
        return true;
    }

    // ===== 按 tstamp 删除（适配未编号对象） =====
    size_t removeViaByTstamp(const std::string& t) {
        size_t n = removeNodesByTstampRec(root, t, "via");
        if (n > 0) rebuildCaches();
        return n;
    }
    size_t removeSegmentByTstamp(const std::string& t) {
        size_t n = removeNodesByTstampRec(root, t, "segment");
        if (n > 0) rebuildCaches();
        return n;
    }

    //存储layer vector
    // 获取顶层板层列表
    const std::vector<LayerInfo>& getBoardLayers() const {
        return boardLayers;
    }

    // 调试打印：板层一览
    void printBoardLayers() const {
        std::cout << "\n=== 板层定义 (共 " << boardLayers.size() << " 层) ===\n";
        for (const auto& L : boardLayers) {
            std::cout << std::setw(2) << L.id << "  "
                << L.name << "  "
                << L.kind;
            if (!L.description.empty()) {
                std::cout << "  \"" << L.description << "\"";
            }
            std::cout << "\n";
        }
    }


    //新增segment

    // 计算下一枚 segment 的编号：现有 [[SEID:x]] 的最大值 + 1；若没有，则从 1 开始
    int nextSegmentId() const {
        int mx = 0;
        for (const auto& s : allSegments) {
            if (!s || s->parameters.empty()) continue;
            int id = extractTagId(s->parameters.front(), "SEID");
            if (id > mx) mx = id;
        }
        return mx + 1;
    }

    // 通用版：layers 可传 1 个（通常 1 个），返回新 segment 的编号（SEID）
    int addSegment(double startX, double startY, double endX, double endY,
        double width, const std::string& layer, int net, const std::string& tstamp)
    {
        if (!root) return -1;

        // 1) 生成新编号与哨兵
        const int newId = nextSegmentId();
        const std::string tag = makeTag("SEID", newId);

        // 2) 构造 segment 节点（结构与原文件一致：segment 自身带哨兵参数，其余都是子节点）
        auto segment = std::make_shared<Node>("segment");
        segment->parameters.push_back(tag);   // [[SEID:x]]

        // (start x y)
        {
            auto startN = std::make_shared<Node>("start");
            startN->parameters.push_back(fmtNum(startX));
            startN->parameters.push_back(fmtNum(startY));
            segment->addChild(startN);
        }

        // (end x y)
        {
            auto endN = std::make_shared<Node>("end");
            endN->parameters.push_back(fmtNum(endX));
            endN->parameters.push_back(fmtNum(endY));
            segment->addChild(endN);
        }

        // (width w)
        {
            auto widthN = std::make_shared<Node>("width");
            widthN->parameters.push_back(fmtNum(width));
            segment->addChild(widthN);
        }

        // (layer)
        {
            auto layerN = std::make_shared<Node>("layer");
            layerN->parameters.push_back(layer);
            segment->addChild(layerN);
        }

        // (net n)
        {
            auto netN = std::make_shared<Node>("net");
            netN->parameters.push_back(std::to_string(net));
            segment->addChild(netN);
        }

        // (tstamp xxx)
        {
            auto tN = std::make_shared<Node>("tstamp");
            tN->parameters.push_back(tstamp);
            segment->addChild(tN);
        }

        // 3) 插到顶层：和原生 segment 一样挂在 root 直下
        root->addChild(segment);

        // 4) 刷新缓存（不重编号，保持其它对象的 ID 稳定）
        rebuildCaches();

        return newId;
    }

    // 便捷版：命令行式输入（层(layer)传空串或“-”表示一个层）
    int addSegmentSimple(double startX, double startY, double endX, double endY,
        double width, const std::string& layer,
        int net, const std::string& tstamp)
    {
        return addSegment(startX, startY, endX, endY, width, layer, net, tstamp);
    }



    // 通用版：layers 可传 0/1/N 个（通常 1 或 2 个），isFree= true/false
// 返回新 via 的编号（VIID）
    int addVia(double x, double y,
        double size, double drill,
        const std::vector<std::string>& layers,
        bool isFree,
        int net,
        const std::string& tstamp)
    {
        if (!root) return -1;

        // 1) 生成新编号与哨兵
        const int newId = nextViaId();
        const std::string tag = makeTag("VIID", newId);

        // 2) 构造 via 节点（结构与原文件一致：via 自身带哨兵参数，其余都是子节点）
        auto via = std::make_shared<Node>("via");
        via->parameters.push_back(tag);   // [[VIID:x]]

        // (at x y)
        {
            auto atN = std::make_shared<Node>("at");
            atN->parameters.push_back(fmtNum(x));
            atN->parameters.push_back(fmtNum(y));
            via->addChild(atN);
        }
        // (size s)
        {
            auto sN = std::make_shared<Node>("size");
            sN->parameters.push_back(fmtNum(size));
            via->addChild(sN);
        }
        // (drill d)
        {
            auto dN = std::make_shared<Node>("drill");
            dN->parameters.push_back(fmtNum(drill));
            via->addChild(dN);
        }
        // (layers ...) —— 有些只有 1 个参数，有些 2 个；我们按传入 vector 动态生成
        if (!layers.empty()) {
            auto layN = std::make_shared<Node>("layers");
            for (const auto& L : layers) {
                if (!L.empty()) layN->parameters.push_back(L); // 内部存不带引号
            }
            // 仅当确实有参数时才插入
            if (!layN->parameters.empty()) via->addChild(layN);
        }
        // (free) 可选
        if (isFree) {
            via->addChild(std::make_shared<Node>("free"));  // 无参数
        }
        // (net n)
        {
            auto nN = std::make_shared<Node>("net");
            nN->parameters.push_back(std::to_string(net));
            via->addChild(nN);
        }
        // (tstamp xxx)
        {
            auto tN = std::make_shared<Node>("tstamp");
            tN->parameters.push_back(tstamp);
            via->addChild(tN);
        }

        // 3) 插到顶层：和原生 via 一样挂在 root 直下
        root->addChild(via);

        // 4) 刷新缓存（不重编号，保持其它对象的 ID 稳定）
        rebuildCaches();

        return newId;
    }

    // 便捷版：按你说的“命令行式”参数（layer2 传空串或 "-" 代表只有 1 个 layer；freeFlag 0/1）
    int addViaSimple(double x, double y,
        double size, double drill,
        const std::string& layer1,
        const std::string& layer2_or_dash,
        int freeFlag,
        int net,
        const std::string& tstamp)
    {
        std::vector<std::string> layers;
        if (!layer1.empty()) layers.push_back(layer1);
        if (!layer2_or_dash.empty() && layer2_or_dash != "-") layers.push_back(layer2_or_dash);
        return addVia(x, y, size, drill, layers, freeFlag != 0, net, tstamp);
    }

    void stripAllTempIds() { stripAllTempIdsRec(root); }

    bool parseFile(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "无法打开文件: " << filename << std::endl;
            return false;
        }

        content.assign((std::istreambuf_iterator<char>(file)),
            std::istreambuf_iterator<char>());
        file.close();

        pos = 0;
        root = parseNode();

        if (root) {
            allNodes.clear();                 // 先清空
            collectAllNodes(root);            // 一次性收集全部节点

            parseBoardLayers();

            findAndNumberSegments();
            findAndNumberVias();
            findAndNumberFootprints();
        }

        return root != nullptr;
    }

    // 获取所有segment
    const std::vector<std::shared_ptr<Node>>& getAllSegments() const {
        return allSegments;
    }

    // 获取所有via
    const std::vector<std::shared_ptr<Node>>& getAllVias() const {
        return allVias;
    }

    // 获取所有footprint
    const std::vector<std::shared_ptr<Node>>& getAllFootprints() const {
        return allFootprints;
    }

    //获取所有节点Node
    const std::vector<std::shared_ptr<Node>>& getAllNodes() const {
        return allNodes;
    }

    // === 对外删除 API ===
// 根据编号删除所有匹配的 via（基于 [[VIID:x]]），返回删除数量
    size_t removeViaById(int id) {
        // 注意：删除依赖哨兵存在，请不要在调用前先 stripAllTempIds()
        size_t n = removeNodesByTagRec(root, "VIID", id, "via");
        if (n > 0) rebuildCaches();
        return n;
    }

    // 根据编号删除所有匹配的 segment（基于 [[SEID:x]]），返回删除数量
    size_t removeSegmentById(int id) {
        size_t n = removeNodesByTagRec(root, "SEID", id, "segment");
        if (n > 0) rebuildCaches();
        return n;
    }


    // 打印segment信息
    void printSegmentInfo(const std::shared_ptr<Node>& segment) {
        if (!segment || segment->name != "segment") {
            std::cout << "不是有效的segment节点" << std::endl;
            return;
        }

        std::cout << "Segment #" << segment->parameters[0] << ":" << std::endl;

        if (segment->parameters.size() > 1) {
            std::cout << "  参数: ";
            for (size_t i = 1; i < segment->parameters.size(); i++) {
                std::cout << "\"" << segment->parameters[i] << "\" ";
            }
            std::cout << std::endl;
        }

        if (!segment->children.empty()) {
            std::cout << "  子节点:" << std::endl;
            for (const auto& child : segment->children) {
                std::cout << "    - " << child->name;
                if (!child->parameters.empty()) {
                    std::cout << " : ";
                    for (const auto& param : child->parameters) {
                        std::cout << "\"" << param << "\" ";
                    }
                }
                std::cout << std::endl;
            }
        }
    }

    // 打印via信息
    void printViaInfo(const std::shared_ptr<Node>& via) {
        if (!via || via->name != "via") {
            std::cout << "不是有效的via节点" << std::endl;
            return;
        }

        // 标题行：显示编号（首参数）
        std::cout << "Via";
        if (!via->parameters.empty()) {
            std::cout << " #" << via->parameters[0];
        }
        std::cout << ":" << std::endl;

        // 参数行：从第二个参数开始打印（跳过编号）
        if (via->parameters.size() > 1) {
            std::cout << "  参数: ";
            for (size_t i = 1; i < via->parameters.size(); ++i) {
                std::cout << "\"" << via->parameters[i] << "\" ";
            }
            std::cout << std::endl;
        }

        // 子节点打印保持不变
        if (!via->children.empty()) {
            std::cout << "  子节点:" << std::endl;
            for (const auto& child : via->children) {
                std::cout << "    - " << child->name;
                if (!child->parameters.empty()) {
                    std::cout << " : ";
                    for (const auto& param : child->parameters) {
                        std::cout << "\"" << param << "\" ";
                    }
                }
                std::cout << std::endl;
            }
        }
    }


    // 打印footprint信息
    void printFootprintInfo(const std::shared_ptr<Node>& footprint) {
        if (!footprint || footprint->name != "footprint") {
            std::cout << "不是有效的footprint节点" << std::endl;
            return;
        }

        std::cout << "Footprint: ";
        if (!footprint->parameters.empty()) {
            std::cout << footprint->parameters[0];
        }
        std::cout << std::endl;

        if (!footprint->children.empty()) {
            std::cout << "  子节点:" << std::endl;
            for (const auto& child : footprint->children) {
                std::cout << "    - " << child->name;
                if (!child->parameters.empty()) {
                    std::cout << " : ";
                    for (const auto& param : child->parameters) {
                        std::cout << "\"" << param << "\" ";
                    }
                }
                std::cout << " (子节点数量: " << child->children.size() << ")" << std::endl;
            }
        }
    }

    // 打印所有segments
    void printAllSegments() {
        if (allSegments.empty()) {
            std::cout << "没有找到任何segment" << std::endl;
            return;
        }

        std::cout << "\n=== 所有Segment列表 (共 " << allSegments.size() << " 个) ===" << std::endl;
        for (const auto& segment : allSegments) {
            printSegmentInfo(segment);
            std::cout << std::endl;
        }
    }

    // 打印所有vias
    void printAllVias() {
        if (allVias.empty()) {
            std::cout << "没有找到任何via" << std::endl;
            return;
        }

        std::cout << "\n=== 所有Via列表 (共 " << allVias.size() << " 个) ===" << std::endl;
        for (const auto& via : allVias) {
            printViaInfo(via);
            std::cout << std::endl;
        }
    }

    // 打印所有footprints
    void printAllFootprints() {
        if (allFootprints.empty()) {
            std::cout << "没有找到任何footprint" << std::endl;
            return;
        }

        std::cout << "\n=== 所有Footprint列表 (共 " << allFootprints.size() << " 个) ===" << std::endl;
        for (const auto& footprint : allFootprints) {
            printFootprintInfo(footprint);
            std::cout << std::endl;
        }
    }

    // 打印树状结构
    void printStructure(const std::shared_ptr<Node>& node = nullptr, int depth = 0, int maxDepth = 3) {
        std::shared_ptr<Node> current = node ? node : root;

        if (!current || depth > maxDepth) {
            return;
        }

        std::string indent(depth * 2, ' ');
        std::cout << indent << "Node: " << current->name;

        if (!current->parameters.empty()) {
            std::cout << " | Parameters: ";
            for (const auto& param : current->parameters) {
                std::cout << "\"" << param << "\" ";
            }
        }
        std::cout << std::endl;

        for (const auto& child : current->children) {
            printStructure(child, depth + 1, maxDepth);
        }
    }
};




// 假设之前的 Node 和 KiCadParser 类都已包含在内，以下是新的主函数

// 递归打印节点结构
void printNodeStructure(const std::shared_ptr<Node>& node, int depth = 0) {
    if (!node) return;



    // 输出当前节点的名称和参数
    std::string indent(depth * 2, ' ');  // 缩进
    std::cout << indent << "节点名称: " << node->name << std::endl;

    if (!node->parameters.empty()) {
        std::cout << indent << "  参数: ";
        for (const auto& param : node->parameters) {
            std::cout << "\"" << param << "\" ";
        }
        std::cout << std::endl;
    }

    // 递归输出子节点
    if (!node->children.empty()) {
        std::cout << indent << "  子节点:" << std::endl;
        for (const auto& child : node->children) {
            printNodeStructure(child, depth + 1);  // 深度加1进行递归
        }
    }
}

// === 辅助：通用的递归打印（用于打印 via 的子树） ===
static void printNodeGeneric(const std::shared_ptr<Node>& node, int depth = 0) {
    if (!node) return;
    std::string indent(depth * 2, ' ');

    // 行首：节点名
    std::cout << indent << "- " << node->name;

    // 同行：节点参数（如果有）
    if (!node->parameters.empty()) {
        std::cout << " : ";
        for (const auto& p : node->parameters) {
            std::cout << "\"" << p << "\" ";
        }
    }
    std::cout << std::endl;

    // 递归打印子节点
    for (const auto& child : node->children) {
        printNodeGeneric(child, depth + 1);
    }
}

// === 专用于 via 的完整递归打印（首参视为编号，与 segment/footprint 一致） ===
static void printViaFull(const std::shared_ptr<Node>& via, int depth = 0) {
    if (!via || via->name != "via") return;
    std::string indent(depth * 2, ' ');

    // 标题：Via #编号
    std::cout << indent << "Via";
    if (!via->parameters.empty()) {
        std::cout << " #" << via->parameters[0]; // 约定：首参为编号
    }
    std::cout << std::endl;

    // 打印除编号外的其它参数（若有）
    if (via->parameters.size() > 1) {
        std::cout << indent << "  参数: ";
        for (size_t i = 1; i < via->parameters.size(); ++i) {
            std::cout << "\"" << via->parameters[i] << "\" ";
        }
        std::cout << std::endl;
    }

    // 递归打印整个子树（包含 at/size/drill/layers/net 等）
    if (!via->children.empty()) {
        std::cout << indent << "  子树:" << std::endl;
        for (const auto& child : via->children) {
            printNodeGeneric(child, depth + 2);
        }
    }
}

#include <iostream>
#include <string>

// 假设 First_Part7.cpp 中已经定义了 KiCadParser 并实现：
// parseFile(...) / removeViaById(...) / removeViaByTstamp(...)
// saveAsKicadPcb(...)

int main() {
    

  
     std::string input  = R"(D:\梁彦诗的项目\Liang的科研\7.0.0版本\C50.kicad_pcb)";
     std::string output = R"(D:\梁彦诗的项目\Liang的科研\混杂的版本\C50.kicad_pcb_after_delete.kicad_pcb)";

    // === 2) 解析文件 ===
    KiCadParser parser;
    if (!parser.parseFile(input)) {
        std::cerr << "[错误] 解析失败： " << input << std::endl;
        return 1;
    }

    // === 3) 删除编号为 3 的 via（按 VIID）===
    // 注意：此操作依赖你的文件已完成 via 编号并带有 [[VIID:3]] 哨兵信息
    size_t removed = parser.removeViaById(3);
    if (removed == 0) {
        std::cerr << "[提示] 未找到 VIID==3 的 via。可能此文件未进行via编号，"
            "或该编号不存在。若只知道 tstamp，可改用：\n"
            "  parser.removeViaByTstamp(\"<该via的tstamp>\");\n";
        // 示例（如果你知道 tstamp 就解除注释使用）：
        // removed = parser.removeViaByTstamp("c8f740da-ca8d-4b56-a5ba-eebb9d2716ee");
    }
    else {
        std::cout << "[信息] 已删除 via (VIID==3) 数量: " << removed << std::endl;
    }

    // === 4) 导出 .kicad_pcb ===
    // 第2个参数 stripTempIds=true：导出前自动清除 [[FPID]]/[[SEID]]/[[VIID]] 哨兵
    if (!parser.saveAsKicadPcb(output, /*stripTempIds=*/true, /*indentStep=*/2)) {
        std::cerr << "[错误] 导出失败： " << output << std::endl;
        return 2;
    }

    std::cout << "[完成] 已生成：" << output << std::endl;
    return 0;
}

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
        if (currentChar() != '"') {
            return "";
        }

        pos++;
        std::string result;

        while (pos < content.size() && currentChar() != '"') {
            if (currentChar() == '\\') {
                pos++;
            }
            if (pos < content.size()) {
                result += currentChar();
                pos++;
            }
        }

        if (currentChar() == '"') {
            pos++;
        }

        return result;
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

        pos++;
        skipWhitespace();

        std::string name;
        while (pos < content.size()) {
            char c = currentChar();
            if (std::isspace(c) || c == '(' || c == ')') {
                break;
            }
            name += c;
            pos++;
        }

        if (name.empty()) {
            return nullptr;
        }

        auto node = std::make_shared<Node>(name);

        while (pos < content.size()) {
            skipWhitespace();

            if (currentChar() == ')') {
                pos++;
                break;
            }

            if (currentChar() == '(') {
                auto child = parseNode();
                if (child) {
                    node->addChild(child);
                }
            }
            /*
            else {
                std::string param;
                if (currentChar() == '"') {
                    param = readQuotedString();
                }
                else {
                    param = readToken();
                }

                if (!param.empty() && param != ")" && param != "(") {
                    node->addParameter(param);
                }
            }
            */
            else {
                std::string param;
                if (currentChar() == '"') {
                    param = readQuotedString();      // 引号内的仍然用它
                }
                else {
                    param = readBareToken();         // 非引号参数用“裸 token”读取
                }
                if (!param.empty()) {
                    node->addParameter(param);
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




//footprint编号调试
/*
int main() {
    std::string inputFile = R"(D:\梁彦诗的项目\Liang的科研\7.0.0版本\A30.kicad_pcb)";

    try {
        // 创建 KiCadParser 实例并解析文件
        KiCadParser parser;
        if (!parser.parseFile(inputFile)) {
            std::cerr << "文件解析失败" << std::endl;
            return 1;
        }

        // 不需要手动调用 findAndNumberFootprints()，因为它在 parseFile 内部已经调用了

        // 获取所有footprint并显示编号
        const auto& allFootprints = parser.getAllFootprints();
        std::cout << "找到 " << allFootprints.size() << " 个footprints" << std::endl;

        // 显示footprint编号和名称
        std::cout << "\n=== Footprint编号列表 ===" << std::endl;
        for (const auto& footprint : allFootprints) {
            // 第一个参数是编号，第二个参数是名称
            if (footprint->parameters.size() >= 2) {
                std::cout << "编号: " << footprint->parameters[0]
                    << " | 名称: " << footprint->parameters[1] << std::endl;
            }
            else if (footprint->parameters.size() >= 1) {
                std::cout << "编号: " << footprint->parameters[0]
                    << " | 名称: (无名)" << std::endl;
            }
            else {
                std::cout << "Footprint没有参数" << std::endl;
            }
        }

    }
    catch (const std::exception& e) {
        std::cerr << "发生异常: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "发生未知异常" << std::endl;
        return 1;
    }

    return 0;
}

*/



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

#include <climits>

#include <iostream>
#include <string>
#include <climits>

// 如果你的头文件/类定义已在同一文件上方，可不需要额外 include
// 
// #include "First_Part4.hpp"  // 仅示例

int main(int argc, char* argv[]) {
    // 默认测试文件路径（按你的实际路径改）
    std::string inputFile = R"(D:\梁彦诗的项目\Liang的科研\混杂的版本\7.0测试版本1.txt)";
    if (argc >= 2) inputFile = argv[1];

    try {
        KiCadParser parser;

        // 1) 解析
        if (!parser.parseFile(inputFile)) {
            std::cerr << "文件解析失败\n";
            return 1;
        }

        // 解析后初始统计
        std::cout << "Before -> vias: " << parser.getAllVias().size()
            << ", segments: " << parser.getAllSegments().size()
            << ", footprints: " << parser.getAllFootprints().size()
            << std::endl;

        // 2) 删除：via #2 与 segment #2
        const int viaIdToRemove = 2;
        const int segIdToRemove = 2;

        size_t removedVia = parser.removeViaById(viaIdToRemove);
        size_t removedSeg = parser.removeSegmentById(segIdToRemove);

        std::cout << "Removed " << removedVia << " via(s) with ID " << viaIdToRemove << "\n";
        std::cout << "Removed " << removedSeg << " segment(s) with ID " << segIdToRemove << "\n";

        // 删除后统计
        std::cout << "After  -> vias: " << parser.getAllVias().size()
            << ", segments: " << parser.getAllSegments().size()
            << ", footprints: " << parser.getAllFootprints().size()
            << std::endl;

        // 3) 打印：保留编号哨兵，按层级输出“allnode”的结构视图
        std::cout << "\n=========== 当前整棵树（保留编号哨兵） ===========" << std::endl;

        // 如果你有 KiCadParser::printStructure(root, depth, maxDepth)
        parser.printStructure(/*root*/nullptr, /*depth*/0, /*maxDepth*/INT_MAX);

        // —— 若你没有 printStructure(...)，而是有全局函数 printNodeStructure(Node*)
        // const auto& allTop = parser.getTopLevelNodes(); // 如果你有这个
        // for (const auto& n : allTop) printNodeStructure(n);

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "发生异常: " << e.what() << "\n";
        return 1;
    }
    catch (...) {
        std::cerr << "发生未知异常\n";
        return 1;
    }
}

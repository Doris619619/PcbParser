#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <memory>
#include <cctype>
#include <algorithm>
#include <functional>


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


    void skipWhitespace() {
        while (pos < content.size() && std::isspace(content[pos])) {
            pos++;
        }
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



    //给所有segment编号
    void findAndNumberSegments() {
        segmentCounter = 0;
        allSegments.clear();
        findNodesRecursive(root, "segment", allSegments);

        for (auto& segment : allSegments) {
            segmentCounter++;
            segment->parameters.insert(segment->parameters.begin(), std::to_string(segmentCounter));
        }
    }



    void findAndNumberVias() {
        viaCounter = 0;
        allVias.clear();
        findNodesRecursive(root, "via", allVias);
        for (auto& via : allVias) {
            ++viaCounter;
            // 把编号插到参数列表最前面（与 segment/footprint 的做法一致）
            via->parameters.insert(via->parameters.begin(), std::to_string(viaCounter));
        }
    }


    void findAndNumberFootprints() {
        int footprintCounter = 0;  // 用来给footprint编号
        allFootprints.clear();  // 清空原有的footprints
        findNodesRecursive(root, "footprint", allFootprints);  // 查找所有footprint节点

        for (auto& footprint : allFootprints) {
            footprintCounter++;  // 递增footprint编号
            footprint->parameters.insert(footprint->parameters.begin(), std::to_string(footprintCounter));  // 插入编号作为第一个参数
        }
    }


public:
    KiCadParser() : pos(0), segmentCounter(0), viaCounter(0) {}

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
            // 查找并存储所有节点
            findNodesRecursive(root, "", allNodes);  // 查找所有节点并存储

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

// === 调试 via 的 main：打印所有 via 的完整信息 ===
int main() {
    std::string inputFile = R"(D:\梁彦诗的项目\Liang的科研\7.0.0版本\A30.kicad_pcb)"; // 修改成你的路径

    try {
        KiCadParser parser;
        if (!parser.parseFile(inputFile)) {
            std::cerr << "文件解析失败" << std::endl;
            return 1;
        }

        const auto& allVias = parser.getAllVias();
        std::cout << "=== 所有 Via (共 " << allVias.size() << " 个) ===" << std::endl;

        for (const auto& via : allVias) {
            printViaFull(via);
            std::cout << std::endl;
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

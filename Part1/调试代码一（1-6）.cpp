/*调试一：调试计算绝对坐标的函数
int main() {
    std::string inputFile = R"(D:\梁彦诗的项目\Liang的科研\7.0.0版本\A30.kicad_pcb)";

    try {
        // 调用计算函数
        std::vector<FootprintPadAbsolute> result = calculateFootprintCoordinates(inputFile);

        std::cout << "解析成功！" << std::endl;
        std::cout << "找到 " << result.size() << " 个pad绝对坐标信息" << std::endl << std::endl;

        // 设置输出格式
        std::cout << std::fixed << std::setprecision(6);

        // 打印所有footprint pad的绝对坐标
        for (const auto& fpa : result) {
            std::cout << "Footprint: " << fpa.footprintName << std::endl;
            std::cout << "  Origin: (" << fpa.footprintOriginX << ", " << fpa.footprintOriginY
                << "), Angle: " << fpa.footprintAngle << "°" << std::endl;
            std::cout << "  Pad: " << fpa.padName
                << " -> Absolute: (" << fpa.padAbsoluteX << ", " << fpa.padAbsoluteY << ")" << std::endl;
            std::cout << "----------------------------------------" << std::endl;
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




//调试二：footprint编号调试
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

//调试三：测试AllNode是否存储所有子节点
//调试三：测试AllNode是否存储所有子节点
// 你已有：KiCadParser / Node / printNodeStructure(...)
/*
int main(int argc, char* argv[]) {
    // ① 文件路径：支持命令行参数覆盖；否则用你的默认路径
    std::string inputFile = R"(D:\梁彦诗的项目\Liang的科研\混杂的版本\7.0测试版本1.txt)";
    if (argc > 1) inputFile = argv[1];

    try {
        KiCadParser parser;
        if (!parser.parseFile(inputFile)) {
            std::cerr << "文件解析失败\n";
            return 1;
        }

        const auto& allNodes = parser.getAllNodes();
        std::cout << "=== 总节点数: " << allNodes.size() << " ===\n";

        if (allNodes.empty()) {
            std::cout << "没有任何节点可打印。\n";
            return 0;
        }

        // ② 计算所有“子节点指针集合”，据此筛出顶层根节点
        std::unordered_set<const void*> childPtrSet;
        childPtrSet.reserve(allNodes.size() * 2);

        for (const auto& n : allNodes) {
            if (!n) continue;
            for (const auto& ch : n->children) {
                if (ch) childPtrSet.insert(static_cast<const void*>(ch.get()));
            }
        }

        std::vector<std::shared_ptr<Node>> roots;
        roots.reserve(4);
        for (const auto& n : allNodes) {
            if (!n) continue;
            if (childPtrSet.find(static_cast<const void*>(n.get())) == childPtrSet.end()) {
                roots.push_back(n); // 只要它从未作为任何人的 child 出现过，就是顶层根
            }
        }

        std::cout << "=== 顶层根节点个数: " << roots.size() << " ===\n\n";

        // ③ 逐棵树打印完整结构（每棵树只打印一次）
        for (size_t i = 0; i < roots.size(); ++i) {
            std::cout << "----- 根 " << (i + 1) << " -----\n";
            printNodeStructure(roots[i]);   // 你已有的递归打印函数
            std::cout << "\n";
        }

        // ④ 附：allNodes 的扁平清单（不递归，便于快速核对）
        std::cout << "=== 扁平清单（名称 | 参数数 | 子节点数）===\n";
        for (const auto& n : allNodes) {
            if (!n) continue;
            std::cout << n->name << " | "
                << n->parameters.size() << " | "
                << n->children.size() << "\n";
        }

    }
    catch (const std::exception& e) {
        std::cerr << "发生异常: " << e.what() << "\n";
        return 1;
    }
    catch (...) {
        std::cerr << "发生未知异常\n";
        return 1;
    }

    return 0;
}
*/


//调试四：打印所有有编号的via
/*
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
*/



//调试五：打印有编号和无编号的allnode
//
/*
// ===================== 放在 First_Part.cpp 文件最下面 =====================
#include <climits>  // 为了使用 INT_MAX

int main(int argc, char* argv[]) {
    // 1) 输入文件：支持命令行参数覆盖
    std::string inputFile = R"(D:\梁彦诗的项目\Liang的科研\混杂的版本\7.0测试版本1.txt)";
    if (argc > 1) inputFile = argv[1];

    try {
        KiCadParser parser;
        if (!parser.parseFile(inputFile)) {
            std::cerr << "文件解析失败\n";
            return 1;
        }

        // A. 先打印——带编号（含哨兵）的“整棵树”，保证层级关系
        std::cout << "\n================ 带编号（含哨兵） ================\n";
        // 传一个极大的 maxDepth，完整 DFS 打印整棵树（从 root 开始）
        parser.printStructure(nullptr, /*indent=*/0, /*maxDepth=*/INT_MAX);

        // B. 去掉哨兵，再打印一次“无编号”的整棵树（同样保留层级）
        parser.stripAllTempIds();  // 全树剥离 [[FPID:]]/[[SEID:]]/[[VIID:]]
        std::cout << "\n================ 无编号（已剥离） ================\n";
        parser.printStructure(nullptr, /*indent=*/0, /*maxDepth=*/INT_MAX);

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



*/


//调试六：删除编号为2的via以及编号为segment的via，并且print出allnode
/*
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

*/

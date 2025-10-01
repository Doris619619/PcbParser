//调试七：新增via
/*
int main() {
    // 输入文件路径（根据你的实际路径调整）
    std::string inputFile = R"(D:\梁彦诗的项目\Liang的科研\混杂的版本\7.0测试版本1.txt)";

    try {
        // 创建 KiCadParser 实例
        KiCadParser parser;

        // 1) 解析文件
        if (!parser.parseFile(inputFile)) {
            std::cerr << "文件解析失败" << std::endl;
            return 1;
        }

        // 2) 打印原始 allnode（无编号，保留结构）
        std::cout << "=== 原始 allnode ===" << std::endl;
        parser.printStructure(nullptr, 0, INT_MAX);  // 打印原始结构（无编号）

        // 3) 新增 via
        int newViaId = parser.addViaSimple(
            149.226102, 115.928602,  // x, y
            0.4, 0.2,                // size, drill
            "F.Cu", "B.Cu",          // layers
            0,                       // freeFlag = 0 (不生成 free 节点)
            35,                      // net
            "c6dc6c1d-fc37-4502-a5ff-45dcc1a3c949" // tstamp
        );

        std::cout << "Added new via with VIID = " << newViaId << std::endl;

        // 4) 再次打印更新后的 allnode（无编号，保留结构）
        std::cout << "=== 更新后的 allnode ===" << std::endl;
        parser.printStructure(nullptr, 0, INT_MAX);  // 打印更新后的结构（无编号）

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









//调试八：新增segment
//我想要添加这个segment： (segment (start 123.140575 105.498365) (end 123.140575 106.628362) (width 1) (layer "B.Cu") (net 39) (tstamp e743a546-2f4c-4132-bb1e-e261cf4cd2f5))
/*
int main(int argc, char* argv[]) {
    // 输入文件路径（根据你的实际路径调整）
    std::string inputFile = R"(D:\梁彦诗的项目\Liang的科研\混杂的版本\7.0测试版本1.txt)";

    try {
        KiCadParser parser;

        // 1) 解析文件
        if (!parser.parseFile(inputFile)) {
            std::cerr << "文件解析失败" << std::endl;
            return 1;
        }

        // 打印原始 allnode（无编号，保留结构）
        std::cout << "=== 原始 allnode ===" << std::endl;
        parser.printStructure(nullptr, 0, INT_MAX);  // 打印原始结构（无编号）

        // 2) 新增 segment
        int newSegmentId = parser.addSegmentSimple(
            123.140575, 105.498365,   // start x, y
            123.140575, 106.628362,   // end x, y
            1,                        // width
            "B.Cu",                   // layer
            39,                       // net
            "e743a546-2f4c-4132-bb1e-e261cf4cd2f5" // tstamp
        );

        std::cout << "Added new segment with SEID = " << newSegmentId << std::endl;

        // 3) 再次打印更新后的 allnode（无编号，保留结构）
        std::cout << "=== 更新后的 allnode ===" << std::endl;
        parser.printStructure(nullptr, 0, INT_MAX);  // 打印更新后的结构（无编号）

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




//调试九：print layer的信息（测试存储layer的vector）
/*
int main(int argc, char* argv[]) {
    // 1) 读取输入文件：优先用命令行参数，否则使用你之前的测试文件
    std::string inputFile;
    if (argc >= 2) {
        inputFile = argv[1];
    }
    else {
        inputFile = R"(D:\梁彦诗的项目\Liang的科研\混杂的版本\7.0测试版本1.txt)";
        // 如果你把文件放在程序当前目录，也可改为：
        // inputFile = "7.0测试版本1.txt";
    }

    try {
        KiCadParser parser;

        // 2) 解析文件
        if (!parser.parseFile(inputFile)) {
            std::cerr << "解析失败：无法读取或解析文件： " << inputFile << std::endl;
            return 1;
        }

        // 3) 只打印板层信息（依赖你刚才添加的 parseBoardLayers / printBoardLayers）
        parser.printBoardLayers();

        // 如需自行使用向量处理，可参考（此段仅示例，保持注释即可）：
        // const auto& layers = parser.getBoardLayers();
        // for (const auto& L : layers) {
        //     // L.id, L.name, L.kind, L.description
        // }

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "发生异常: " << e.what() << std::endl;
        return 2;
    }
}

*/



//调试十：输出kicad文件（并且删除编号为3的via）（示例里面用的是C50那个kicad文件作为输入）
/*
int main() {
    

  
     std::string input  = R"(D:\梁彦诗的项目\Liang的科研\混杂的版本\7.0测试版本1.txt)";//改为你的地址
     std::string output = R"(D:\梁彦诗的项目\Liang的科研\混杂的版本\7.0测试版本1_after_delete.kicad_pcb)";

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

*/

//调试板框信息，输出板框长宽
//输入文件是MIJI-LCM-EL922-926.kicad_pcb
/*
int main(int argc, char** argv) {
    // 1) 从命令行取文件路径；2) 没给参数就用本地测试文件，自己改成 L92.txt/A30.txt/B40.txt 的路径
    std::string inputFile;
    if (argc >= 2) {
        inputFile = argv[1];
    }
    else {
        inputFile = R"(D:\梁彦诗的项目\Liang的科研\7.0.0版本\L92.txt)"; // ← 改成你的样例路径
    }

    try {
        KiCadParser parser;
        if (!parser.parseFile(inputFile)) {
            std::cerr << "解析失败: " << inputFile << std::endl;
            return 1;
        }

        // 计算外接矩形与长宽
        KiCadParser::BoardBBox bbox = parser.getBoardBBox();
        if (!bbox.valid) {
            std::cout << "未在 Edge.Cuts 找到板框（gr_line/gr_arc/gr_rect/gr_poly）。" << std::endl;
            return 0;
        }

        double width = 0.0, height = 0.0;
        (void)parser.getBoardSize(width, height); // 与 bbox.width()/height() 等价

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "文件: " << inputFile << "\n";
        std::cout << "板框外接矩形 (单位: mm)\n";
        std::cout << "  minX = " << bbox.minX << ", minY = " << bbox.minY << "\n";
        std::cout << "  maxX = " << bbox.maxX << ", maxY = " << bbox.maxY << "\n";
        std::cout << "  宽度(width)  = " << width << "\n";
        std::cout << "  高度(height) = " << height << "\n";

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "异常: " << e.what() << std::endl;
        return 2;
    }
}

*/

//调试十二：获取net信息（类比layer）
/*
    int main() {
    std::string inputFile = R"(D:\梁彦诗的项目\Liang的科研\7.0.0版本\L92.txt)";

    try {
        KiCadParser parser;

        if (!parser.parseFile(inputFile)) {
            std::cerr << "解析失败：无法读取文件或语法错误 -> " << inputFile << std::endl;
            return 1;
        }

        const auto& nets = parser.getBoardNets();

        std::cout << "=== 网络表 (共 " << nets.size() << " 条) ===\n";
        for (const auto& n : nets) {
            std::cout << "(net " << n.id << " \"" << n.name << "\")\n";
        }

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "异常：" << e.what() << std::endl;
        return 2;
    }
}
*/

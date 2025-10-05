调试一：输入layer 0，区域编号1，输出该阻塞区域内的所有阻塞像素的像素编号（左下角第一个像素是（1,1））

  // ---------------- 示例 main ----------------
// 用法：./your_app [kicad_text_file]
// 不传参时默认读取 "A30.txt"（你也可以改成 L92.txt / testcase.txt 等）
int main(int argc, char** argv) {
    try {
        std::string filename = (argc >= 2) ? argv[1] : "testcase.kicad_pcb";

        Grid grid;
        grid.SetUp(filename);  // 解析 KiCad 文本并构建网格、离散化等（你们已有流程）:contentReference[oaicite:1]{index=1}

        int layer = 0;
        int areaID = 1;

        PixelList1 cells1 = GetBlockedAreaPixels_1Based(grid, layer, areaID);

        std::cout << "Layer " << layer << ", AreaID " << areaID
            << " 的像素数量 = " << cells1.size() << std::endl;

        // 打印所有像素（1-based）
        for (const auto& rc1 : cells1) {
            std::cout << "(row=" << rc1.first << ", col=" << rc1.second << ")";
            // 同时演示：该像素中心对应的现实坐标（mm）
            auto p = PixelCenterToReal_1Based(grid, rc1.first, rc1.second);
            std::cout << " -> center(x=" << p.first << ", y=" << p.second << ")\n";
        }

        // 演示：现实坐标 → 像素(1-based)
        {
            double sx = 100.02183, sy = 91.97056;
            auto rc1 = RealToPixel_1Based(grid, sx, sy);
            std::cout << "示例坐标 (" << sx << ", " << sy
                << ") 所在像素(1-based) = (row=" << rc1.first
                << ", col=" << rc1.second << ")\n";
        }
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "[Divide.cpp] 运行出错: " << e.what() << std::endl;
        return 1;
    }
}

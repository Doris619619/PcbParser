#include "grid.h"
#include <cmath>
#include <algorithm>
#include <queue>

// ��ȫ��ȡһ�� segment/via/pad �� net �ַ����������Ҳ����ͷ��ؿմ���
static std::string GetNetIdOfNode(const Node* n) {
    if (!n) return "";
    for (const auto& ch : n->children) {
        if (ch && ch->name == "net" && !ch->parameters.empty()) {
            return ch->parameters[0];        // ע�⣺�������ַ���ID��"5"������������
        }
    }
    return "";
}


void Grid::SetUp(std::string filename) {
    KiCadParser parser;
    // ʹ�� KiCadParser �����ļ�

    if (parser.parseFile(filename)) {
        std::cout << "�����ɹ���" << std::endl;
    }
    else {
        std::cerr << "����ʧ�ܣ�" << std::endl;
        return;
    }
    footprints = parser.getAllFootprints();  // ȷ�����ݴ��ݸ� grid
    vias = parser.getAllVias();  // ��ȡ���� vias
    segments = parser.getAllSegments();
    parser.printAllSegments();
    parser.printAllFootprints();
    std::vector<NetInfo> nets = parser.getBoardNets();
    KiCadParser::Point2D bl, tl, tr, br;
    if (!parser.getBoardCorners(bl, tl, tr, br)) {
        std::cout << "[δ�ҵ� Edge.Cuts ���] " << std::endl;
        return;
    }
    for (auto net : nets) {
        NetNamewithID.insert({ net.id,net.name });
        //std::cout<<net.id<<" "<<net.name;
    }

    // ��ȡ���� segments
    for (auto layerstruct : parser.getBoardLayers()) {
        layerName.push_back(layerstruct.name);
    }
    absCoor = calculateFootprintCoordinates(filename);
    MatchLayers();
    find_max_min_coor(bl, tr);
    initialize();
    fill_grid();
    CheckSpecialBlock();
    PartitionRegions();
}
int Grid::getLayerId(const std::string& layerName) {
    if (layerIDwithName.find(layerName) != layerIDwithName.end()) {
        return layerIDwithName[layerName];  // ���ز�� ID
    }
    else {
        std::cerr << "Error: Layer name not found: " << layerName << std::endl;
        return -1;  // ����Ҳ������������ش���ֵ
    }
}

void Grid::initialize() {
    Layers = layerIDwithName.size();
    width = int((max_x - min_x) * inputScale);
    height = int((max_y - min_y) * inputScale);
    grid.resize(Layers, Grid2D(height, std::vector<cell>(width)));
}

// ���������С����
void Grid::find_max_min_coor(KiCadParser::Point2D bl, KiCadParser::Point2D tr) {
    min_x = bl.x;
    min_y = bl.y;
    max_x = tr.x;
    max_y = tr.y;
    //std::cout<<min_x<<' '<<min_y<<' '<<max_x<<' '<<max_y;
}

void Grid::SetMinSize(int minSize) {
    minSize = minSize;
}
void Grid::fill_grid() {
    for (auto& fp : footprints) {
        for (auto child : fp->children) {
            if (child->name == "pad") {
                //double pad_grid_x = FindPadAbsoluteCoor(stod(fp->children[2]->parameters[0]), stod(fp->children[2]->parameters[1]), child->parameters[0])[0] - min_x;
                //double pad_grid_y = FindPadAbsoluteCoor(stod(fp->children[2]->parameters[0]), stod(fp->children[2]->parameters[1]), child->parameters[0])[1] - min_y;


                //�޸�
                auto getFpAt = [&](const std::shared_ptr<Node>& fp, double& x, double& y)->bool {
                    for (const auto& ch : fp->children) {
                        if (ch && ch->name == "at" && ch->parameters.size() >= 2) {
                            x = std::stod(ch->parameters[0]);
                            y = std::stod(ch->parameters[1]);
                            return true;
                        }
                    }
                    return false;
                    };

                double fpX = 0, fpY = 0;
                if (!getFpAt(fp, fpX, fpY)) {
                    // �Ҳ��� at����������� footprint �� pad������Խ��
                    continue;
                }





                // ��ȡ footprint ԭ����Ƕ�
                double fpAng = 0.0;
                for (const auto& ch : fp->children) {
                    if (ch && ch->name == "at" && ch->parameters.size() >= 2) {
                        fpX = std::stod(ch->parameters[0]);
                        fpY = std::stod(ch->parameters[1]);
                        if (ch->parameters.size() >= 3) fpAng = std::stod(ch->parameters[2]);
                        break;
                    }
                }

                // ���� Calculator �Ľ��
                auto padAbs = FindPadAbsoluteCoor(fpX, fpY, child->parameters[0]);
                if (padAbs.size() < 2) {
                    std::cerr << "[Warn] pad abs not found by Calculator, fallback to Node tree\n";
                    // ���� ���ˣ�ֱ���ýڵ�����һ�Σ�������B��ļ��㹫ʽ��
                    double dx = 0, dy = 0;
                    for (const auto& z : child->children) {
                        if (z->name == "at" && z->parameters.size() >= 2) {
                            dx = std::stod(z->parameters[0]);
                            dy = std::stod(z->parameters[1]);
                            break;
                        }
                    }
                    double th = fpAng * M_PI / 180.0;
                    double padXAbs = fpX + std::cos(th) * dx - std::sin(th) * dy;
                    double padYAbs = fpY + std::sin(th) * dx + std::cos(th) * dy;
                    padAbs = { padXAbs, padYAbs };
                }

                double pad_grid_x = padAbs[0] - min_x;
                double pad_grid_y = padAbs[1] - min_y;


                //double pad_grid_x = FindPadAbsoluteCoor(fpX, fpY, child->parameters[0])[0] - min_x;
                //double pad_grid_y = FindPadAbsoluteCoor(fpX, fpY, child->parameters[0])[1] - min_y;

                float radius;
                float width;
                float height;
                int angle = 0;
                std::string NetName;
                std::vector<int> layersID;
                bool isGND = false;
                for (auto Child : child->children) {
                    if (Child->name == "layer" || Child->name == "layers") {
                        // �滻Ϊ�µ�ͨ����������߼�
                        auto addLayerToken = [&](const std::string& token) {
                            auto addIf = [&](const std::string& name) {
                                auto it = layerIDwithName.find(name);
                                if (it != layerIDwithName.end()) layersID.push_back(it->second);
                                };

                            if (token == "*.Cu" || token == "F&B.Cu") {          // ͨ��ͭ�� / ˫��ͭ
                                addIf("F.Cu"); addIf("B.Cu");
                            }
                            else if (token == "*.Mask") {                       // ͨ���躸
                                addIf("F.Mask"); addIf("B.Mask");
                            }
                            else if (token == "*.Paste") {                      // ͨ������
                                addIf("F.Paste"); addIf("B.Paste");
                            }
                            else if (token == "*.SilkS") {                      // ͨ��˿ӡ
                                addIf("F.SilkS"); addIf("B.SilkS");
                            }
                            else {
                                int id = getLayerId(token);                       // ��ͨ����
                                if (id >= 0 && id < Layers) layersID.push_back(id);
                            }
                            };

                        for (const std::string& token : Child->parameters) {
                            addLayerToken(token);
                        }

                        // ȥ�أ�����ͬһ���ظ�����
                        std::sort(layersID.begin(), layersID.end());
                        layersID.erase(std::unique(layersID.begin(), layersID.end()), layersID.end());
                    }
                    else if (Child->name == "size") {
                        width = std::stod(Child->parameters[0]);
                        height = std::stod(Child->parameters[1]);
                        radius = std::stod(Child->parameters[0]) / 2;
                    }
                    else if (Child->name == "at") {
                        if (Child->parameters.size() > 2) {
                            angle = std::stoi(Child->parameters[2]);
                        }
                    }
                    else if (Child->name == "net") {
                        NetName = NetNamewithID.at(stoi(Child->parameters[0]));
                        if (NetName.find("GND") != std::string::npos) {
                            isGND = true;
                        }
                    }

                }
                std::string pad_shape = child->parameters[child->parameters.size() - 1];
                if (pad_shape == "rect") {
                    for (int layerId : layersID) {
                        if (layerId < 0 || layerId >= Layers) continue; // ����
                        if (isGND) {
                            rectToGrid(Grid::width, Grid::height, pad_grid_x, pad_grid_y, width, height, angle, layerId, inputScale, grid, &(*child));
                        }
                        else {
                            rectToGrid(Grid::width, Grid::height, pad_grid_x, pad_grid_y, width + Clearance, height + Clearance, angle, layerId, inputScale, grid, &(*child));
                        }

                    }
                }
                else if (pad_shape == "circle") {

                    for (int layerId : layersID) {
                        if (layerId < 0 || layerId >= Layers) continue; // ����
                        if (isGND) {
                            circleToGrid(Grid::width, Grid::height, pad_grid_x, pad_grid_y, radius, layerId, inputScale, grid, &(*child));
                        }
                        else {
                            circleToGrid(Grid::width, Grid::height, pad_grid_x, pad_grid_y, radius + Clearance, layerId, inputScale, grid, &(*child));
                        }
                    }
                }
                else if (pad_shape == "roundrect") {
                    float rratio = 0.0;
                    for (auto Child : child->children) {
                        if (Child->name == "roundrect_rratio") {
                            // �޸ģ�ȥ���ظ��� float ����
                            rratio = 0.0f;
                            for (const auto& ch : child->children) {
                                if (ch && ch->name == "roundrect_rratio" && !ch->parameters.empty()) {
                                    rratio = std::stof(ch->parameters[0]);
                                    break;
                                }
                            }
                        }
                    }

                    for (int layerId : layersID) {
                        if (layerId < 0 || layerId >= Layers) continue; // ����
                        if (isGND) {
                            RoundRectToGrid(Grid::width, Grid::height, pad_grid_x, pad_grid_y, width, height, rratio, angle, layerId, inputScale, grid, &(*child));
                        }
                        else {
                            RoundRectToGrid(Grid::width, Grid::height, pad_grid_x, pad_grid_y, width + Clearance, height + Clearance, rratio, angle, layerId, inputScale, grid, &(*child));
                        }
                    }

                }
                else if (pad_shape == "oval") {

                    for (int layerId : layersID) {
                        if (layerId < 0 || layerId >= Layers) continue; // ����
                        if (isGND) {
                            ovalToGrid(Grid::width, Grid::height, pad_grid_x, pad_grid_y, width, height, angle, layerId, inputScale, grid, &(*child));
                        }
                        else {
                            ovalToGrid(Grid::width, Grid::height, pad_grid_x, pad_grid_y, width + Clearance, height + Clearance, angle, layerId, inputScale, grid, &(*child));
                        }
                    }

                }
            }
        }
    }

    for (auto& via : vias) {
        double via_grid_x;
        double via_grid_y;
        float via_radius;
        std::vector<int> layersID;
        for (std::shared_ptr<Node>& child : via->children) {
            if (child->name == "layers" || child->name == "layer") {
                for (auto& param : child->parameters) {
                    layersID.push_back(layerIDwithName.at(param));
                }
            }
            else if (child->name == "at") {
                via_grid_x = std::stod(child->parameters[0]) - min_x;
                via_grid_y = std::stod(child->parameters[1]) - min_y;
            }
            else if (child->name == "size") {
                via_radius = std::stod(child->parameters[0]) / 2;
            }
        }
        for (int layerId : layersID) {
            if (layerId < 0 || layerId >= Layers) continue; // ����
            circleToGrid(Grid::width, Grid::height, via_grid_x, via_grid_y, via_radius, layerId, inputScale, grid, &(*via));
        }
    }

    for (auto& segment : segments) {
        double x1, y1, x2, y2;
        float width;
        std::string NetName;
        std::vector<int> layersID;
        bool isGND = false;
        for (std::shared_ptr<Node>& child : segment->children) {
            if (child->name == "layer") {
                for (auto& param : child->parameters) {
                    layersID.push_back(layerIDwithName.at(param));
                }
            }
            else if (child->name == "start") {
                x1 = std::stod(child->parameters[0]) - min_x;
                y1 = std::stod(child->parameters[1]) - min_y;
            }
            else if (child->name == "end") {
                x2 = std::stod(child->parameters[0]) - min_x;
                y2 = std::stod(child->parameters[1]) - min_y;
            }
            else if (child->name == "width") {
                width = std::stof(child->parameters[0]);
            }
            else if (child->name == "net") {
                NetName = NetNamewithID.at(stoi(child->parameters[0]));
                if (NetName.find("GND") != std::string::npos) {
                    isGND = true;
                    //std::cout<<NetName<<' ';
                }
            }
        }
        for (auto& layerId : layersID) {
            if (layerId < 0 || layerId >= Layers) continue; // ����
            if (isGND) {
                lineToGrid(isGND, Grid::width, Grid::height, x1, y1, x2, y2, width, 0, layerId, inputScale, grid, &(*segment));
            }
            else {
                lineToGrid(isGND, Grid::width, Grid::height, x1, y1, x2, y2, width, Clearance, layerId, inputScale, grid, &(*segment));
            }
        }
    }


}



// ���������������
void Grid::setInputScale(int w) {
    inputScale = w;
}

// �������ֵ
void Grid::setClearance(float clearance) {
    Clearance = clearance;
}

// ��ƥ��
void Grid::MatchLayers() {
    // = {"F.Cu", "B.Cu", "B.Adhes", "F.Adhes", "B.Paste", "F.Paste", "B.SilkS", "F.SilkS", "B.Mask", "F.Mask", "Dwgs.User", "Cmts.User", "Eco1.User", "Eco2.User", "Edge.Cuts", "Margin", "B.CrtYd", "F.CrtYd", "B.Fab", "F.Fab", "User.1", "User.2", "User.3", "User.4", "User.5", "User.6", "User.7", "User.8", "User.9","Rescue"};
    for (int i = 0; i < layerName.size(); i++) {
        layerIDwithName.insert({ layerName[i], i });
    }
}

void Grid::CheckSpecialBlock() {
    for (int i = 0; i < Layers; i++) {
        for (int j = 0; j < height; j++) {
            for (int k = 0; k < width; k++) {
                if (grid[i][j][k].isblocked == false) {
                    for (auto node1 : grid[i][j][k].Segment_Components) {
                        if (k - 1 >= 0 && grid[i][j][k - 1].isblocked == false) {
                            for (auto node2 : grid[i][j][k - 1].Segment_Components) {
                                if (GetNetIdOfNode(node1) != GetNetIdOfNode(node2)) {
                                    grid[i][j][k].isblocked = true;
                                    grid[i][j][k - 1].isblocked = true;
                                    grid[i][j][k].Segment_Components.push_back(node2);
                                    grid[i][j][k - 1].Segment_Components.push_back(node1);
                                    break;
                                }
                            }
                            if (grid[i][j][k].isblocked == true) {
                                break;
                            }
                        }
                        if (k + 1 <= width - 1 && grid[i][j][k + 1].isblocked == false) {
                            for (auto node2 : grid[i][j][k + 1].Segment_Components) {
                                if (GetNetIdOfNode(node1) != GetNetIdOfNode(node2)) {
                                    grid[i][j][k].isblocked = true;
                                    grid[i][j][k + 1].isblocked = true;
                                    grid[i][j][k].Segment_Components.push_back(node2);
                                    grid[i][j][k + 1].Segment_Components.push_back(node1);
                                    break;
                                }
                            }
                            if (grid[i][j][k].isblocked == true) {
                                break;
                            }
                        }
                        if (j - 1 >= 0 && grid[i][j - 1][k].isblocked == false) {
                            for (auto node2 : grid[i][j - 1][k].Segment_Components) {
                                if (GetNetIdOfNode(node1) != GetNetIdOfNode(node2)) {
                                    grid[i][j][k].isblocked = true;
                                    grid[i][j - 1][k].isblocked = true;
                                    grid[i][j][k].Segment_Components.push_back(node2);
                                    grid[i][j - 1][k].Segment_Components.push_back(node1);
                                    break;
                                }
                            }
                            if (grid[i][j][k].isblocked == true) {
                                break;
                            }
                        }
                        if (j + 1 <= height - 1 && grid[i][j + 1][k].isblocked == false) {
                            for (auto node2 : grid[i][j + 1][k].Segment_Components) {
                                if (GetNetIdOfNode(node1) != GetNetIdOfNode(node2)) {
                                    grid[i][j][k].isblocked = true;
                                    grid[i][j + 1][k].isblocked = true;
                                    grid[i][j][k].Segment_Components.push_back(node2);
                                    grid[i][j + 1][k].Segment_Components.push_back(node1);
                                    break;
                                }
                            }
                            if (grid[i][j][k].isblocked == true) {
                                break;
                            }
                        }
                        if (k - 1 >= 0 && j - 1 >= 0 && grid[i][j - 1][k - 1].isblocked == false) {
                            for (auto node2 : grid[i][j - 1][k - 1].Segment_Components) {
                                if (GetNetIdOfNode(node1) != GetNetIdOfNode(node2)) {
                                    grid[i][j][k].isblocked = true;
                                    grid[i][j - 1][k - 1].isblocked = true;
                                    grid[i][j][k].Segment_Components.push_back(node2);
                                    grid[i][j - 1][k - 1].Segment_Components.push_back(node1);
                                    break;
                                }
                            }
                            if (grid[i][j][k].isblocked == true) {
                                break;
                            }
                        }
                        if (k - 1 >= 0 && j + 1 <= height - 1 && grid[i][j + 1][k - 1].isblocked == false) {
                            for (auto node2 : grid[i][j + 1][k - 1].Segment_Components) {
                                if (GetNetIdOfNode(node1) != GetNetIdOfNode(node2)) {
                                    grid[i][j][k].isblocked = true;
                                    grid[i][j + 1][k - 1].isblocked = true;
                                    grid[i][j][k].Segment_Components.push_back(node2);
                                    grid[i][j + 1][k - 1].Segment_Components.push_back(node1);
                                    break;
                                }
                            }
                            if (grid[i][j][k].isblocked == true) {
                                break;
                            }
                        }
                        if (k + 1 <= width - 1 && j - 1 >= 0 && grid[i][j - 1][k + 1].isblocked == false) {
                            for (auto node2 : grid[i][j - 1][k + 1].Segment_Components) {
                                if (GetNetIdOfNode(node1) != GetNetIdOfNode(node2)) {
                                    grid[i][j][k].isblocked = true;
                                    grid[i][j - 1][k + 1].isblocked = true;
                                    grid[i][j][k].Segment_Components.push_back(node2);
                                    grid[i][j - 1][k + 1].Segment_Components.push_back(node1);
                                    break;
                                }
                            }
                            if (grid[i][j][k].isblocked == true) {
                                break;
                            }
                        }
                        if (k + 1 <= width - 1 && j + 1 <= height - 1 && grid[i][j + 1][k + 1].isblocked == false) {
                            for (auto node2 : grid[i][j + 1][k + 1].Segment_Components) {
                                if (GetNetIdOfNode(node1) != GetNetIdOfNode(node2)) {
                                    grid[i][j][k].isblocked = true;
                                    grid[i][j + 1][k + 1].isblocked = true;
                                    grid[i][j][k].Segment_Components.push_back(node2);
                                    grid[i][j + 1][k + 1].Segment_Components.push_back(node1);
                                    break;
                                }
                            }
                            if (grid[i][j][k].isblocked == true) {
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}

/*
std::vector<double> Grid::FindPadAbsoluteCoor(double fp_x, double fp_y, std::string padname) {
    std::vector<double> coor;
    for (auto abscoor : absCoor) {
        if (abscoor.footprintOriginX == fp_x && abscoor.footprintOriginY == fp_y && abscoor.padName == padname) {
            coor.push_back(abscoor.padAbsoluteX);
            coor.push_back(abscoor.padAbsoluteY);
        }
    }
    return coor;
}
*/

// grid.cpp
std::vector<double> Grid::FindPadAbsoluteCoor(double fp_x, double fp_y, std::string padname) {
    std::vector<double> coor;
    const double EPS = 1e-6;
    for (const auto& abscoor : absCoor) {
        if (std::fabs(abscoor.footprintOriginX - fp_x) < EPS &&
            std::fabs(abscoor.footprintOriginY - fp_y) < EPS &&
            abscoor.padName == padname) {
            coor.push_back(abscoor.padAbsoluteX);
            coor.push_back(abscoor.padAbsoluteY);
            break; // �ҵ��͹���
        }
    }
    return coor;
}


void rectToGrid(int max_x, int max_y, double centerX, double centerY, float width, float height, int angle, int layer, int inputScale, std::vector<std::vector<std::vector<cell>>>& grid, Node* temp) {
    if (angle == 90 || angle == -90 || angle == 270) {
        rectToGrid(max_x, max_y, centerX, centerY, height, width, 0, layer, inputScale, grid, temp);
    }
    else {
        for (int j = std::max(int((centerY - height / 2) * inputScale), 0); j <= std::min(int((centerY + height / 2) * inputScale), max_y); j++) {
            for (int i = std::max(int((centerX - width / 2) * inputScale), 0); i <= std::min(int((centerX + width / 2) * inputScale), max_x); i++) {
                bool check = false;
                for (auto node : grid[layer][j][i].nodes) {
                    if (node == temp) {
                        check = true;
                    }
                }
                if (check == true) {
                    continue;
                }
                grid[layer][j][i].isoccupied = 1;
                grid[layer][j][i].nodes.push_back(temp);

            }
        }
    }
}


void circleToGrid(int max_x, int max_y, double centerX, double centerY, float radius, int layer, int inputScale, std::vector<std::vector<std::vector<cell>>>& grid, Node* temp) {
    for (int j = std::max(int((centerY - radius) * inputScale), 0); j <= std::min(int((centerY + radius) * inputScale), max_y); j++) {
        for (int i = std::max(int((centerX - radius) * inputScale), 0); i <= std::min(int((centerX + radius) * inputScale), max_x); i++) {
            const double nx = std::min(std::max(centerX * inputScale, (double)i), (double)i + 1.0) - centerX * inputScale;
            const double ny = std::min(std::max(centerY * inputScale, (double)j), (double)j + 1.0) - centerY * inputScale;
            if (nx * nx + ny * ny <= (radius * inputScale) * (radius * inputScale)) {
                bool check = false;
                for (auto node : grid[layer][j][i].nodes) {
                    if (node == temp) {
                        check = true;
                    }
                }
                if (check == true) {
                    continue;
                }
                grid[layer][j][i].isoccupied = 1;
                grid[layer][j][i].nodes.push_back(temp);


            }
        }
    }
}


void lineToGrid(bool isGND, int max_x, int max_y, double x1, double y1, double x2, double y2, float width, float clearance, int layer, int inputScale, std::vector<std::vector<std::vector<cell>>>& grid, Node* temp) {
    if (x1 == x2) {
        for (int j = std::max(int(std::min(y1, y2) * inputScale), 0); j <= std::min(int(std::max(y1, y2) * inputScale), max_y); j++) {
            for (int i = std::max(int((x1 - width / 2 - clearance) * inputScale), 0); i <= std::min(int((x1 + width / 2 + clearance) * inputScale), max_x); i++) {
                if (isGND) {
                    grid[layer][j][i].isoccupied = 1;
                    grid[layer][j][i].nodes.push_back(temp);
                    continue;
                }
                grid[layer][j][i].nodes.push_back(temp);
                grid[layer][j][i].Segment_Components.push_back(temp);
                if (j == 0 || j == max_y || i == 0 || i == max_x) {
                    grid[layer][j][i].isblocked = 1;
                }
                if (grid[layer][j][i].isoccupied == 1) {
                    for (auto node : grid[layer][j][i].Segment_Components) {
                        if (GetNetIdOfNode(node) != GetNetIdOfNode(temp)) {
                            grid[layer][j][i].isblocked = 1;
                            break;
                        }
                    }
                }
                else {
                    grid[layer][j][i].isoccupied = 1;

                }
            }
        }
    }
    else if (y1 == y2) {
        for (int j = std::max(int((y1 - width / 2 - clearance) * inputScale), 0); j <= std::min(int((y1 + width / 2 + clearance) * inputScale), max_y); j++) {
            for (int i = std::max(int(std::min(x1, x2) * inputScale), 0); i <= std::min(int(std::max(x1, x2) * inputScale), max_x); i++) {
                if (isGND) {
                    grid[layer][j][i].isoccupied = 1;
                    grid[layer][j][i].nodes.push_back(temp);
                    continue;
                }
                grid[layer][j][i].nodes.push_back(temp);
                grid[layer][j][i].Segment_Components.push_back(temp);
                if (j == 0 || j == max_y || i == 0 || i == max_x) {
                    grid[layer][j][i].isblocked = 1;
                }
                if (grid[layer][j][i].isoccupied == 1) {
                    for (auto node : grid[layer][j][i].Segment_Components) {
                        if (GetNetIdOfNode(node) != GetNetIdOfNode(temp)) {
                            grid[layer][j][i].isblocked = 1;
                            break;
                        }
                    }
                }
                else {
                    grid[layer][j][i].isoccupied = 1;

                }
            }
        }
    }
    else {
        double minX = std::min(x1, x2) - clearance;
        double maxX = std::max(x1, x2) + clearance;
        double minY = std::min(y1, y2) - clearance;
        double maxY = std::max(y1, y2) + clearance;

        //45��
        if ((y2 - y1) / (x2 - x1) > 0) {
            for (int j = std::max(int(minY * inputScale), 0); j <= std::min(int(maxY * inputScale), max_y); j++) {
                for (int i = std::max(int(minX * inputScale), 0); i <= std::min(int(maxX * inputScale), max_x); i++) {
                    double t1 = ((i / inputScale - x1) * (x2 - x1) + (j / inputScale - y1) * (y2 - y1)) / ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
                    double t4 = (((i + 1) / inputScale - x1) * (x2 - x1) + ((j + 1) / inputScale - y1) * (y2 - y1)) / ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
                    if ((t1 > 1 || t1 < 0) && (t4 > 1 || t4 < 0)) {
                        continue;
                    }

                    if ((j + 1 - y1 * inputScale) - (y2 - y1) / (x2 - x1) * (i - x1 * inputScale) > 0 && (j - y1 * inputScale) - (y2 - y1) / (x2 - x1) * (i + 1 - x1 * inputScale) < 0) {
                        if (isGND) {
                            grid[layer][j][i].isoccupied = 1;
                            grid[layer][j][i].nodes.push_back(temp);
                            continue;
                        }
                        grid[layer][j][i].nodes.push_back(temp);
                        grid[layer][j][i].Segment_Components.push_back(temp);
                        if (j == 0 || j == max_y || i == 0 || i == max_x) {
                            grid[layer][j][i].isblocked = 1;
                        }
                        if (grid[layer][j][i].isoccupied == 1) {
                            for (auto node : grid[layer][j][i].Segment_Components) {
                                if (GetNetIdOfNode(node) != GetNetIdOfNode(temp)) {
                                    grid[layer][j][i].isblocked = 1;
                                    break;
                                }
                            }
                        }
                        else {
                            grid[layer][j][i].isoccupied = 1;

                        }
                        continue;
                    }

                    double A = y2 - y1;
                    double B = x1 - x2;
                    double C = x2 * y1 - x1 * y2;
                    double distanceTL = std::abs(A * i / inputScale + B * (j + 1) / inputScale + C) / std::sqrt(A * A + B * B);
                    double distanceBR = std::abs(A * (i + 1) / inputScale + B * (j + 1) / inputScale + C) / std::sqrt(A * A + B * B);

                    if (distanceTL <= width / 2 + clearance || distanceBR <= width / 2 + clearance) {
                        if (isGND) {
                            grid[layer][j][i].isoccupied = 1;
                            grid[layer][j][i].nodes.push_back(temp);
                            continue;
                        }
                        grid[layer][j][i].nodes.push_back(temp);
                        grid[layer][j][i].Segment_Components.push_back(temp);
                        if (j == 0 || j == max_y || i == 0 || i == max_x) {
                            grid[layer][j][i].isblocked = 1;
                        }
                        if (grid[layer][j][i].isoccupied == 1) {
                            for (auto node : grid[layer][j][i].Segment_Components) {
                                if (GetNetIdOfNode(node) != GetNetIdOfNode(temp)) {
                                    grid[layer][j][i].isblocked = 1;
                                    break;
                                }
                            }
                        }
                        else {
                            grid[layer][j][i].isoccupied = 1;

                        }
                    }
                }
            }
        }
        else if ((y2 - y1) / (x2 - x1) < 0) {
            for (int j = std::max(int(minY * inputScale), 0); j <= std::min(int(maxY * inputScale), max_y); j++) {
                for (int i = std::max(int(minX * inputScale), 0); i <= std::min(int(maxX * inputScale), max_x); i++) {
                    double t2 = (((i + 1) / inputScale - x1) * (x2 - x1) + (j / inputScale - y1) * (y2 - y1)) / ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
                    double t3 = ((i / inputScale - x1) * (x2 - x1) + ((j + 1) / inputScale - y1) * (y2 - y1)) / ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
                    if ((t2 > 1 || t2 < 0) && (t3 > 1 || t3 < 0)) {
                        continue;
                    }
                    if ((j + 1 - y1 * inputScale) - (y2 - y1) / (x2 - x1) * (i + 1 - x1 * inputScale) > 0 && (j - y1 * inputScale) - (y2 - y1) / (x2 - x1) * (i - x1 * inputScale) < 0) {
                        if (isGND) {
                            grid[layer][j][i].isoccupied = 1;
                            grid[layer][j][i].nodes.push_back(temp);
                            continue;
                        }
                        grid[layer][j][i].nodes.push_back(temp);
                        grid[layer][j][i].Segment_Components.push_back(temp);
                        if (j == 0 || j == max_y || i == 0 || i == max_x) {
                            grid[layer][j][i].isblocked = 1;
                        }
                        if (grid[layer][j][i].isoccupied == 1) {
                            for (auto node : grid[layer][j][i].Segment_Components) {
                                if (GetNetIdOfNode(node) != GetNetIdOfNode(temp)) {
                                    grid[layer][j][i].isblocked = 1;
                                    break;
                                }
                            }
                        }
                        else {
                            grid[layer][j][i].isoccupied = 1;

                        }
                        continue;
                    }
                    double A = y2 - y1;
                    double B = x1 - x2;
                    double C = x2 * y1 - x1 * y2;
                    double distanceTR = std::abs(A * (i + 1) / inputScale + B * (j + 1) / inputScale + C) / std::sqrt(A * A + B * B);
                    double distanceBL = std::abs(A * i / inputScale + B * j / inputScale + C) / std::sqrt(A * A + B * B);
                    if (distanceTR <= width / 2 + clearance || distanceBL <= width / 2 + clearance) {
                        if (isGND) {
                            grid[layer][j][i].isoccupied = 1;
                            grid[layer][j][i].nodes.push_back(temp);
                            continue;
                        }
                        grid[layer][j][i].nodes.push_back(temp);
                        grid[layer][j][i].Segment_Components.push_back(temp);
                        if (j == 0 || j == max_y || i == 0 || i == max_x) {
                            grid[layer][j][i].isblocked = 1;
                        }
                        if (grid[layer][j][i].isoccupied == 1) {
                            for (auto node : grid[layer][j][i].Segment_Components) {
                                if (GetNetIdOfNode(node) != GetNetIdOfNode(temp)) {
                                    grid[layer][j][i].isblocked = 1;
                                    break;
                                }
                            }
                        }
                        else {
                            grid[layer][j][i].isoccupied = 1;

                        }
                    }
                }
            }
        }
    }
    if (isGND) {
        LineCircleToGrid(isGND, max_x, max_y, x1, y1, width / 2, layer, inputScale, grid, temp);
        LineCircleToGrid(isGND, max_x, max_y, x2, y2, width / 2, layer, inputScale, grid, temp);
    }
    else {
        LineCircleToGrid(isGND, max_x, max_y, x1, y1, width / 2 + clearance, layer, inputScale, grid, temp);
        LineCircleToGrid(isGND, max_x, max_y, x2, y2, width / 2 + clearance, layer, inputScale, grid, temp);
    }
}

void LineCircleToGrid(bool isGND, int max_x, int max_y, double centerX, double centerY, float radius, int layer, int inputScale, std::vector<std::vector<std::vector<cell>>>& grid, Node* temp) {
    for (int j = std::max(int((centerY - radius) * inputScale), 0); j <= std::min(int((centerY + radius) * inputScale), max_y); j++) {
        for (int i = std::max(int((centerX - radius) * inputScale), 0); i <= std::min(int((centerX + radius) * inputScale), max_x); i++) {
            const double nx = std::min(std::max(centerX * inputScale, (double)i), (double)i + 1.0) - centerX * inputScale;
            const double ny = std::min(std::max(centerY * inputScale, (double)j), (double)j + 1.0) - centerY * inputScale;
            if (nx * nx + ny * ny <= (radius * inputScale) * (radius * inputScale)) {
                bool check = false;
                for (auto node : grid[layer][j][i].nodes) {
                    if (node == temp) {
                        check = true;
                    }
                }
                if (check == true) {
                    continue;
                }
                if (isGND) {
                    grid[layer][j][i].isoccupied = 1;
                    grid[layer][j][i].nodes.push_back(temp);
                    continue;
                }
                grid[layer][j][i].nodes.push_back(temp);
                grid[layer][j][i].Segment_Components.push_back(temp);
                if (j == 0 || j == max_y || i == 0 || i == max_x) {
                    grid[layer][j][i].isblocked = 1;
                }
                if (grid[layer][j][i].isoccupied == 1) {
                    for (auto node : grid[layer][j][i].Segment_Components) {
                        if (GetNetIdOfNode(node) != GetNetIdOfNode(temp)) {
                            grid[layer][j][i].isblocked = 1;
                            break;
                        }
                    }
                }
                else {
                    grid[layer][j][i].isoccupied = 1;

                }

            }
        }
    }
}
void RoundRectToGrid(int max_x, int max_y, double CenterX, double CenterY, float Width, float Height, float rratio, int Angle, int Layer, int InputScale, std::vector<std::vector<std::vector<cell>>>& Grid, Node* Temp) {
    float radius = std::min(Width, Height) * rratio;
    //if (chamfer_ratio == 0){
    if (Angle == 90 || Angle == -90 || Angle == 270) {
        RoundRectToGrid(max_x, max_y, CenterX, CenterY, Height, Width, rratio, 0, Layer, InputScale, Grid, Temp);
    }
    else {
        rectToGrid(max_x, max_y, CenterX, CenterY, Width - radius * 2, Height, 0, Layer, InputScale, Grid, Temp);
        rectToGrid(max_x, max_y, CenterX, CenterY, Width, Height - radius * 2, 0, Layer, InputScale, Grid, Temp);
        circleToGrid(max_x, max_y, CenterX - Width / 2 + radius, CenterY - Height / 2 + radius, radius, Layer, InputScale, Grid, Temp);
        circleToGrid(max_x, max_y, CenterX + Width / 2 - radius, CenterY - Height / 2 + radius, radius, Layer, InputScale, Grid, Temp);
        circleToGrid(max_x, max_y, CenterX - Width / 2 + radius, CenterY + Height / 2 - radius, radius, Layer, InputScale, Grid, Temp);
        circleToGrid(max_x, max_y, CenterX + Width / 2 - radius, CenterY + Height / 2 - radius, radius, Layer, InputScale, Grid, Temp);
    }
}

void ovalToGrid(int max_x, int max_y, double centerX, double centerY, float width, float height, int angle, int layer, int inputScale, std::vector<std::vector<std::vector<cell>>>& grid, Node* temp) {
    if (height > width) {
        if (angle == 90 || angle == -90 || angle == 270) {
            ovalToGrid(max_x, max_y, centerX, centerY, height, width, 0, layer, inputScale, grid, temp);
        }
        else {
            /*grid[layer][int(centerY*inputScale)][int(centerX*inputScale)].isoccupied = 1;
            grid[layer][int(centerY*inputScale)][int(centerX*inputScale)].nodes.push_back(temp);

            for(int j = int((centerY - (height-width)/2)*inputScale); j <= int((centerY + (height-width)/2)*inputScale); j++){
                for(int i = int((centerX - width/2)*inputScale); i <= int((centerX + width/2)*inputScale); i++){
                    grid[layer][j][i].isoccupied = 1;
                    grid[layer][j][i].nodes.push_back(temp);

                }
            }*/
            rectToGrid(max_x, max_y, centerX, centerY, width, height - width, 0, layer, inputScale, grid, temp);
            circleToGrid(max_x, max_y, centerX, centerY - ((height - width) / 2), width / 2, layer, inputScale, grid, temp);
            circleToGrid(max_x, max_y, centerX, centerY + ((height - width) / 2), width / 2, layer, inputScale, grid, temp);
        }
    }
    else {
        if (angle == 90 || angle == -90 || angle == 270) {
            ovalToGrid(max_x, max_y, centerX, centerY, height, width, 0, layer, inputScale, grid, temp);
        }
        else {
            /*grid[layer][int(centerY*inputScale)][int(centerX*inputScale)].isoccupied = 1;
            grid[layer][int(centerY*inputScale)][int(centerX*inputScale)].nodes.push_back(temp);

            for(int i = int((centerY - height/2)*inputScale); i <= int((centerY + height/2)*inputScale); i++){
                for(int j = int((centerX - (width-height)/2)*inputScale); j <= int((centerX + (width-height)/2)*inputScale); j++){
                    grid[layer][i][j].isoccupied = 1;
                    grid[layer][i][j].nodes.push_back(temp);

                }
            }*/
            rectToGrid(max_x, max_y, centerX, centerY, width - height, height, 0, layer, inputScale, grid, temp);
            circleToGrid(max_x, max_y, centerX - (width - height) / 2, centerY, height / 2, layer, inputScale, grid, temp);
            circleToGrid(max_x, max_y, centerX + (width - height) / 2, centerY, height / 2, layer, inputScale, grid, temp);
        }
    }
}

void Grid::PartitionRegions() {
    auto inb = [&](int r, int c) -> bool {
        return r >= 0 && r < height && c >= 0 && c < width;
        };

    const int dr[4] = { -1, 1, 0, 0 }; // ��������
    const int dc[4] = { 0, 0,-1, 1 };

    for (int i = 0; i < Layers; ++i) {
        for (int r = 0; r < height; ++r) {
            for (int c = 0; c < width; ++c) {

                // ֻ��"������δ���"�ĸ������� BFS
                if (grid[i][r][c].isoccupied != 0 || grid[i][r][c].regionID != 0)
                    continue;

                // ����ͨ�飺�ȷ���һ����ʱ���
                ++RegionID;
                std::queue<std::pair<int, int>> q;
                std::vector<std::pair<int, int>> cells; // ��¼�ÿ�������и��ӣ����ڻ���
                int count = 0;

                q.push({ r, c });
                grid[i][r][c].regionID = RegionID;

                // BFS ����
                while (!q.empty()) {
                    auto [x, y] = q.front(); q.pop();
                    cells.push_back({ x, y });
                    ++count;

                    for (int k = 0; k < 4; ++k) {
                        int nx = x + dr[k];
                        int ny = y + dc[k];
                        if (!inb(nx, ny)) continue;
                        if (grid[i][nx][ny].isoccupied != 0) continue;
                        if (grid[i][nx][ny].regionID != 0) continue;
                        grid[i][nx][ny].regionID = RegionID;
                        q.push({ nx, ny });
                    }
                }

                // �ж���С��ģ��С����ֵ�����ϸñ�Ų�����
                if (count < minSize) {
                    for (auto& p : cells) {
                        grid[i][p.first][p.second].regionID = 0;
                    }
                    --RegionID; // ���˱�ż����������¿պ�
                }
                // �������� RegionID���������⴦��
            }
        }
    }
}
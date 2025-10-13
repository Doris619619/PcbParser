#ifndef GRID_H
#define GRID_H

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>
#include <stack>
#include "First_Part12.h"
#include "Calculator.h"
#include <queue>
#include <cstdint>

struct cell {
public:
    int isoccupied = 0;
    bool isblocked = false;
    int regionID = 0;

    int cost = 0;

    std::vector<Node*> nodes;
    std::vector<Node*> Segment_Components;
};

class Grid {
public:
    int width;
    int height;
    int Layers;
    float Clearance = 0.7;
    int RegionID = 0;
    int minSize = 20;
    std::vector<FootprintPadAbsolute> absCoor;
    using Grid2D = std::vector<std::vector<cell>>;
    using Grid3D = std::vector<Grid2D>;
    std::vector<std::shared_ptr<Node>> segments;
    std::vector<std::shared_ptr<Node>> vias;
    std::vector<std::shared_ptr<Node>> footprints;
    std::vector<std::string> layerName;
    std::map<std::string, int> layerIDwithName;
    std::map<int, std::string> NetNamewithID;
    Grid3D grid;
    int inputScale = 10;
    double max_x = std::numeric_limits<double>::max();
    double max_y = std::numeric_limits<double>::max();
    double min_x = std::numeric_limits<double>::min();
    double min_y = std::numeric_limits<double>::min();
    std::vector<double> FindPadAbsoluteCoor(double fp_x, double fp_y, std::string padname);

    // Functions for grid initialization and layer handling
    void initialize();
    void find_max_min_coor(KiCadParser::Point2D bl, KiCadParser::Point2D tr);
    void fill_grid();
    void setInputScale(int w);
    void setClearance(float clearance);
    void MatchLayers();
    void MatchNets();
    int getLayerId(const std::string& layerName);
    void PartitionRegions();
    void CheckSpecialBlock();
    void SetUp(std::string filename);
    void SetMinSize(int minSize);

};

void RoundRectToGrid(int max_x, int max_y, double CenterX, double CenterY, float Width, float Height, float rratio, int Angle, int Layer, int InputScale, std::vector<std::vector<std::vector<cell>>>& Grid, Node* Temp);
void lineToGrid(bool isGND, int max_x, int max_y, double x1, double y1, double x2, double y2, float width, float clearance, int layer, int inputScale, std::vector<std::vector<std::vector<cell>>>& grid, Node* temp);
void ovalToGrid(int max_x, int max_y, double centerX, double centerY, float width, float height, int angle, int layer, int inputScale, std::vector<std::vector<std::vector<cell>>>& grid, Node* temp);
void rectToGrid(int max_x, int max_y, double centerX, double centerY, float width, float height, int angle, int layer, int inputScale, std::vector<std::vector<std::vector<cell>>>& grid, Node* temp);
void circleToGrid(int max_x, int max_y, double centerX, double centerY, float radius, int layer, int inputScale, std::vector<std::vector<std::vector<cell>>>& grid, Node* temp);
void rectToGrid(int max_x, int max_y, double centerX, double centerY, float width, float height, int angle, int layer, int inputScale, std::vector<std::vector<std::vector<cell>>>& grid, Node* temp);
void LineCircleToGrid(bool isGND, int max_x, int max_y, double centerX, double centerY, float radius, int layer, int inputScale, std::vector<std::vector<std::vector<cell>>>& grid, Node* temp);

#endif // GRID_H
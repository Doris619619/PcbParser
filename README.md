# Ground-Plane-Generation

## 1. Convert kicad_pcb to kicad_class  



### 《第一次修改》
#### 【说明】目前数据的结构如下：
每一个节点叫Node，节点的参数叫Parameter，如果一个节点长这样子
  (segment (start 148.613828 87.077563) (end 148.513826 86.977561) (width 0.6) (layer "B.Cu") (net 33) (tstamp 7146d6c6-40cf-415e-880d-b8b906abdfa0))
那么segment就是Node，segment的子节点就是start，end，width......，子节点start的parameter就是148.613828 87.077563

#### 【可实现的功能】目前Fisrt_Part.cpp能实现的功能有：
###### 1，存储：将所有via，segment和footprint分别存储在三个vector中。（Vector名字分别是allSegments，allVias，allFootprints;）
###### 2，获取：可以获取这三个vector。（函数分别是getAllSegments，getAllVias，getAllFootprints）
###### 3，编号：给所有segment编号（按照出现顺序编号）
###### 4，检验：比如可以打印出所有包含via的所有节点信息（从而验证前面存储via的vector的正确性）

#### 【注意事项】
###### 1,文件路径：运行代码的时候需要修改文件路径






### 《第二次修改》
#### 【增加的功能】
###### 1，增加了Calculator.cpp文件和Calculator.h文件用以计算pad的绝对坐标。并且把最终结果用名字叫result的一个容器装了起来，容器里装了很多个结构体对象。（// result 是一个vector实例→std::vector<FootprintPadAbsolute> result;）
###### 2，在First_Part.cpp文件中增加了调试Calculator算出来坐标的批注，启用这个批注可以print出所有pad的绝对坐标。（注意要把以前的main函数批注掉）

### 《第三次修改》
###### 【增加的文件】增加了几个调试代码的结果，以便于清晰看到数据的结构



### 《第四次修改》9.29
#### 【增加的文件】
###### 1，调试后的结果文件（调试x.txt）
###### 2，一个包含实现各种调试的main函数的文件（调试文件.cpp），选择不同main函数可以看到不同的输出
###### 3，更新了First_Part.cpp，最新版本的文件是First_Part3.cpp。这个版本相比于First_Part.cpp，不仅修复了层级之间的关系（以前具有一点漏洞会导致读取存在一些bug），而且还新增了如下功能。
#### 【First_Part3.cpp新增的功能】
###### ①给所有footprint和via编号
###### ②增加了一个获取所有数据（包含解析完文件之后的所有node）const std::vector<std::shared_ptr<Node>>& getAllNodes() const {  return allNodes;}
###### ③给所有编号做上哨兵标记。（这样使得最后在打印所有数据的时候是不带上标记的via segment footprint编号的，但是在删除segment和via的时候可以用编号作为输入来删除）
###### gpt解释：继续把编号写到 parameters[0]，但用“哨兵标记”，最后一键剥离。思路：插入带独特前缀的编号（比如 [[FPID:123]]、[[SEID:5]]、[[VIID:9]]），操作完后全树遍历把这些带前缀的参数删除或还原。
###### ④删除segment和via。对于这两个函数，输入是编号。
gpt解释如下：
删除的使用方法很简单：解析文件后（parser.parseFile(path)），via 和 segment 会各自带上稳定编号哨兵（[[VIID:x]] / [[SEID:x]]）。要删除某个对象，只需调用 parser.removeViaById(id) 或 parser.removeSegmentById(id)；它会在整棵语法树中把编号精确等于该值的对应节点全部移除，并自动重建缓存（allNodes/allVias/allSegments），无需手动刷新。例如删除编号为 2 的 via 和 segment：parser.removeViaById(2); parser.removeSegmentById(2); 然后用你的层级打印（如 parser.printStructure(nullptr,0,INT_MAX)）或遍历 getAllNodes() 验证即可。注意：删除依赖这些哨兵存在，所以在删除前不要调用会清理哨兵的函数（如 stripAllTempIds()）；删除不会重新编号，若需重排编号请另行调用你的重编号逻辑。若采用提供的命令行 main，也可以运行：prog <pcb_file> del via 2 或 prog <pcb_file> del segment 2。

#### 【关于调试文件的说明】
###### （所有调试的输入均用的是某个kicad的截取文件，文件名：7.0测试版本1.txt，可在文件夹里找到）
###### ①调试六：实现了删除整个文件中编号为2的via和segment，并且输出删除via和segment之后的allnode。（如果想要调用删除函数，可以看调试六的main函数）


### 《第五次修改》9.30
#### 【First_Part4.cpp新增的功能】（相比于version 3）
###### ①添加via和segment。
要新增 segment 和 via，只需调用以下方法：

新增 segment：
调用 addSegmentSimple(startX, startY, endX, endY, width, layer, net, tstamp)。
输入参数：起点和终点坐标 (startX, startY, endX, endY)，宽度 (width)，层 (layer)，网络编号 (net)，时间戳 (tstamp)。

新增 via：
调用 addViaSimple(x, y, size, drill, layer1, layer2, freeFlag, net, tstamp)。
输入参数：坐标 (x, y)，尺寸 (size)，钻孔大小 (drill)，层 (layer1, layer2)，是否自由 (freeFlag)，网络编号 (net)，时间戳 (tstamp)。

这两个函数会自动为每个 segment 和 via 分配唯一编号，并将其插入到树结构中，确保在 allNodes 和相关缓存中可以找到新节点。
###### ②一个获取layer信息的函数
要获取 layer 信息，可以通过 KiCadParser 类中的 getBoardLayers() 函数，它返回一个包含所有板层信息的 std::vector<LayerInfo>。
每个 LayerInfo 对象存储以下字段：
id：层的编号（如 0, 31, 44 等）
name：层的名称（如 "F.Cu", "B.Cu"）
kind：层的类型（如 "signal", "user"）
description：层的描述（如 "Top", "Bottom"），可选。
使用时，通过 parser.getBoardLayers() 可以访问所有层的信息。

### 《第六次修改》9.30
#### 【First_Part5.cpp新增的功能】（相比于version 4）
①实现了kicad文件的逆转。

### 《第七次修改》10.1
###### ①添加头文件First_Part6.h(为了配合头文件First_Part.cpp也做出了微调，class里面的每一个函数都加上了KiCadParser::）
将 Node/LayerInfo/KiCadParser 的类型定义与全部函数声明抽离到 First_Part10.h，只保留必要的内联小函数，对外统一接口并避免多翻译单元的重复定义风险。
First_Part6.cpp 仅保留 #include "First_Part6.h" 与所有成员函数的类外实现（使用 KiCadParser::作用域限定，函数体与原注释不变），同时删除原先类内实现及重复的结构体/类定义

### 《第八次修改》10.1
#### 【First_Part7.cpp新增的功能】（相比于version6）
###### ①获取板框
现在可通过 KiCadParser::getBoardSize(width, height) 直接获得板框长宽
###### ②获取net信息
获取网络表：调用 getBoardNets()，返回 const std::vector<NetInfo>&。
NetInfo.id 是 (net <id> "...") 中的编号（int）。
NetInfo.name 是网络名（可能是空字符串，比如 (net 0 "")）。




## 2. Initial region partition  
## 3. Generate GPG graph based on routing graph  
## 4. ILP formulation and solving with Gurobi  
## 5. Rerouting






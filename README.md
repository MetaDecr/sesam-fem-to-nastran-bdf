# fem2bdf — GeniE/Sesam FEM 转 Nastran BDF 转换器

> 该项目由 Kimi K3 生成。

## 这个程序是做什么的

GeniE（DNV 的船体与海洋结构建模工具）可以把有限元模型导出为 Sesam
接口文件（俗称 T 文件，扩展名通常是 `.fem` 或 `.FEM`）。这是一种纯文本、
固定格式的卡片文件，里面对节点、单元、截面、材料、边界条件、载荷各有
一套 Sesam 自己的卡片定义（GNODE、GELMNT1、GIORH、BEUSLO……）。

Nastran 不认这套格式。它要的是 Bulk Data File（`.bdf`），同样按卡片组织，
但卡片体系完全不同（GRID、CBEAM、PBEAML、SPC1、PLOAD4……）。

fem2bdf 做的事情就是把前者翻译成后者：读入 GeniE 导出的 `.fem`，
写出一份可以直接交给 Nastran 求解的 SOL 101（线性静力）完整 deck——
包含 Case Control（SUBCASE / SPC / LOAD）和全部 Bulk Data。
转换结束时会在控制台打印统计信息和所有警告（跳过了什么、近似了什么、
缺了什么引用），让用户对转换质量一目了然。

解析和写出都基于 DNV 在 GitHub 上开源的 C++ 读写库
[ExtFEMIO](https://github.com/dnvgl/ExtFEMIO)（源码已内置于
`extfemio/`，有 4 处必要补丁，见文末）。也就是说，这不是一个自己用
字符串切分硬啃固定格式的脚本，而是建立在 DNV 官方维护的格式解析器
之上，卡片的字段边界、续行、大字段格式都由库来处理。

## 转换原理

### 总体流程

程序先把整个 FEM 文件读进一个内存模型（节点、单元、各类截面卡、材料、
边界、载荷分门别类存放），然后扫描每个单元通过 GELREF1 引用的几何号
（GEONO）和材料号（MATNO），只把**真正被单元用到**的属性写出来。
最后按 Nastran 的习惯顺序输出：MAT1 → GRID → 属性卡（PSHELL / PBEAM /
PBEAML / PSOLID / PBUSH……）→ 单元卡 → 边界与载荷。

属性 ID 的分配规则是：同一 GEONO 下第一种材料沿用 GEONO 作为 PID；
同一 GEONO 搭配不同材料时自动分配新 PID，保证不与其他 ID 冲突。
实体单元的 PSOLID PID 由材料号派生。

### 梁截面：形状优先，通用数据兜底

这是本转换器最重要的设计决策之一。GeniE 导出梁截面时有个习惯：
同一个 GEONO 下**同时**写一张 GBEAMG（预先算好的面积、惯性矩等通用
特性）和一张形状卡（GIORH、GLSEC、GPIPE 等描述截面几何尺寸的卡片）。
如果转换器先看到 GBEAMG 就用它，形状信息就丢了。

fem2bdf 的策略是**形状卡优先**：只要存在形状卡，就转换成对应的
PBEAML 库截面（Nastran 按尺寸自己计算截面特性，还能自动给出应力恢复
点）；GBEAMG→PBEAM 的直接特性拷贝只作为没有形状卡时的回退。
球扁钢、H 型钢等 GeniE 标准库截面在接口文件里只以 GBEAMG 输出，
这类就只能走 PBEAM 路径（面积和惯性矩精确，但没有截面形状纤维信息）。

各形状卡的映射关系：

| Sesam (FEM) | Nastran (BDF) | 说明 |
|---|---|---|
| GBEAMG | PBEAM | A=AREA, I1=IY, I2=IZ, J=IX；IYZ≠0 时给 I12=−IYZ；K1/K2 由 SHARZ/SHARY÷AREA 得到。仅在无形状卡时使用 |
| GPIPE | PBEAML (TUBE) | 圆管；DIM = [外径, 壁厚] |
| GIORH | PBEAML (I) | 工字钢/H 型钢；DIM = [HZ, BB, BT, TY, TB, TT] |
| GIORHR | PBEAML (I) | 带内圆角工字钢；圆角半径不保留并警告 |
| GUSYI | PBEAML (I) | 非对称工字钢；腹板偏心（B1/B2≠半宽）忽略并警告 |
| GBOX | PBEAML (BOX/BOX1) | 箱形梁、方管/矩形管；上下翼缘等厚写 BOX，否则写 BOX1 |
| GDOBO | PBEAML (I) | 双底截面；上下板作等宽翼缘（宽 BY） |
| GCHAN / GCHANR | PBEAML (CHAN) | 槽钢；DIM = [BY, HZ, TZ, TY]；GCHANR 圆角忽略并警告 |
| GLSEC / GLSECR | PBEAML (L) | 角钢；DIM = [BY, HZ, TZ, TY]（BY 水平肢、HZ 竖直肢）；GLSECR 圆角忽略并警告 |
| GTONP | PBEAML (I) | T 型材（带板）；有效板（BP×TP）作下翼缘、面板（BT×TT）作上翼缘，总高取 HZ+TP |
| GBARM | PBEAML (BAR) / PBEAM | 扁钢/方钢；等宽矩形写 BAR；梯形截面无对应 Nastran 形状，回退为等效 PBEAM 并警告 |

PBEAML 的 DIM 顺序已按 MSC Nastran Reference Guide 逐一核对
（例如 I 型是 [高, 下翼宽, 上翼宽, 腹板厚, 下翼厚, 上翼厚]，
L 型是 [水平肢宽, 竖直肢高, 水平肢厚, 竖直肢厚]）。

### 梁的偏心：形心轴与剪心轴的修正

这是最容易被忽略、也最容易做错的一点。

Sesam 的梁单元节点位于**形心轴**（center of gravity line）上，GECCEN
卡片给出的是全局坐标系下"网格节点 → 单元节点（形心）"的偏心向量。
而 Nastran 的 CBEAM + PBEAML 是从**剪心**（shear centre）开始绘制
截面的——单元轴就是剪心轴，形心相对剪心的偏移由 Nastran 根据截面
形状自己计算。

对于双轴对称截面（I、管、箱、矩形）两者重合，无所谓；但对角钢、
槽钢这类非对称截面，形心和剪心相差可能超过截面高度的一半。
如果直接把 Sesam 的偏心向量抄成 Nastran 的 WA/WB，整个截面就会
平移到一个错误的位置。

fem2bdf 的处理是：WA/WB = GECCEN + 形心→剪心向量。（局部坐标系偏心卡
GECC 也支持：先按单元局部轴旋转到全局，再走同一通道。）后者取自
GBEAMG 卡片的 SHCENY/SHCENZ 字段（Sesam 局部 y/z 系下的分量），
再按 GUNIVEC 定义的局部坐标系旋转到全局。以一根挂在板下的
L150×90×9 角钢为例：GeniE 给的偏心是形心深度 99.72 mm，修正后
Nastran 的 WA 指向板下 145.5 mm 的肢背脚跟（两翼中线交点）——
这正是角钢剪心的理论位置。对称截面 SHCEN 为零，不受影响；
即使单元没有 GECCEN，只要截面非对称也会自动补上这个修正。

梁的方向向量取 GUNIVEC：Sesam 局部 z 轴对应 Nastran 局部 y 轴。

### 单元与其他卡片映射

| Sesam (FEM) | Nastran (BDF) | 说明 |
|---|---|---|
| GNODE + GCOORD | GRID | 使用外部节点号（NODEX≠NODENO 时按映射转换） |
| MISOSEL | MAT1 | E、NU、RHO；G 留空由 Nastran 自行推导 |
| GELTH | PSHELL | T = 板厚 |
| FQUS / FFQ (24) | CQUAD4 | |
| FTRS / FFTR (25) | CTRIA3 | |
| BEAS (15) / BEPS (2) | CBEAM | 见上面两节 |
| TEXT / TDSETNAM + GSETMEMB | $ 注释 | 用户文本写入文件头注释；集合（名称、类型、成员，节点集用外部编号）作为注释保留在 ENDDATA 之前 |
| TESS (10) | CROD | PROD 由 GBEAMG 的 AREA 或 GPIPE 截面积得到 |
| LHEX (21) / IHEX (20) | CHEXA + PSOLID | 角点序按 det J 检查并自动翻转；20 节点单元按边端点对映射中边节点 |
| TETR (33) / ITET (31) | CTETRA + PSOLID | 同上 |
| TPRI (32) / IPRI (30) | CPENTA + PSOLID | 同上 |
| GSPR (18) + MGSPRNG | CBUSH（接地）+ PBUSH | 取刚度矩阵对角 6 项；非对角耦合项丢弃并警告 |
| MORSMEL | MAT2 | 各向异性膜材料：D 矩阵直接映射 G11~G33，ALPHA1/2→A1/A2；方向向量 Q 无法表示，丢弃并警告 |
| GMAS (11) + MGMASS | CONM1 | 完整 6×6 对称质量矩阵（MGMASS 按列存下三角，已按文档解包） |
| BNMASS | CMASS2 | 每个非零自由度一张卡，EID 自动顺延 |
| BELFIX | CBEAM 释放标记 PA/PB | A(i)=0（完全释放）→ 对应自由度写入 PA/PB；OPT=SPRING 的弹性铰刚度和 0<A<1 的部分固定无法表示，按固接处理并警告；TRANO 局部系忽略并警告 |
| BLDEP | MPC | 每个相关自由度一张 MPC（u_dep − Σbᵢ·uᵢ = 0），子工况自动加 MPC=1；BNBCD 标记为相关（3）但缺 BLDEP 的自由度会警告 |
| BNBCD (FIX) | SPC1 | 按自由度组合分组 |
| BNBCD (prescribed) + BNDISPL | SPC | 每个载荷工况一套 SPC 集（SID = 工况号），DTYPE=1 取值；DTYPE=3（加速度）按 0 约束并警告 |
| BNLOAD | FORCE / MOMENT | 兼容 GeniE 实际布局与官方文档布局两种写法 |
| BEUSLO (LOTYP=±1) | PLOAD4 | 壳：SIDE=1 正压、SIDE=2/3 反号；实体：SIDE 数字旗标定位面号，面序按外法线几何校验 |
| TDLOAD | SUBCASE LABEL | 载荷工况名称 |

### 缺材料的兜底

GeniE 模型如果没有指定材料，接口文件里就不会有 MISOSEL 卡，
但单元引用还在。直接写出 MID=0 的 BDF 是非法的，Nastran 会拒绝。
fem2bdf 遇到这种情况会自动补一张**全零占位 MAT1**（E=ν=ρ=0）
并明确警告。补零而不是补一组看似合理的默认值（比如钢材），是
刻意的：占位值不可能被遗忘——零刚度的模型一求解就报错，
绝不会悄悄算出一个"看起来没问题"的错误结果。

### 缺厚度 / 缺截面的兜底

类似地，板单元引用的 GEONO 可能根本没有 GELTH 卡（或 TH≤0），
即"模型里有板但没有厚度"。此时若简单跳过 PSHELL，壳单元就会
引用一个不存在的 PID，BDF 直接报废。fem2bdf 的做法是照常写出
PSHELL，厚度补零占位，并打印醒目的 WARNING——与材料一样，
零厚度的板无法求解，提醒你必须补上真实厚度。

梁/杆单元引用的 GEONO 如果没有任何截面卡（GBEAMG/GPIPE/GIORH
等都没有），则连占位都做不到（截面有十几个相互关联的参数，
补零会立刻毁掉求解器的单元检查）：该属性不写出，引用它的单元也
一并跳过并统计警告，保证 BDF 中不会出现悬空的属性引用。

## 局限性

转换器追求的是"忠实 + 显式"：能精确转的就精确转，必须近似或放弃的
一律在日志里警告，绝不留沉默的坑。目前已知的局限如下。

**不转换的内容（跳过并统计警告）：**

- 阻尼类单元（GDAM、GLDA 等）与阻尼矩阵卡片。
- BEUSLO 中 LOTYP=2（分量式面载荷）与 LOTYP=3（水动压力标记）；
  LOTYP=−1（随动压力）按普通法向压力转换并警告。

**近似处理的内容：**

- 弹簧/质量单元若定义在局部坐标系（GELREF1 带 TRANSNO），矩阵按
  全局系直接使用并警告；弹簧刚度矩阵的非对角耦合项丢弃
  （PBUSH 只保留 6 个对角项）。
- 梯形 GBARM 截面回退为等效 PBEAM（面积/惯性矩精确，无应力恢复点）。
- GUSYI 非对称工字钢的腹板偏心忽略，按腹板居中处理并警告。
- GIORHR / GLSECR / GCHANR 的内圆角半径不保留，按直角截面近似并警告。
- GTONP（带板 T 型材）的高度约定：HZ 按 T 部高度处理，总高 = HZ+TP，
  有效板整体作为下翼缘写入 I 型 DIM。
- GDOBO 双底截面按上下板等宽（BY）的 I 型近似。
- BELFIX 铰接只有"完全释放"（A=0）能精确映射为 CBEAM 释放标记；
  部分固定（0<A<1）和 OPT=SPRING 的弹性铰刚度无法表示，按固接处理
  并警告；TRANO 局部系下的铰接按单元局部系处理并警告。
- MORSMEL 各向异性膜材料的方向向量 Q 无法表示（MAT2 不携带方向），
  丢弃并警告；DAMP1/2、PS1/2 也不转换。
- BLDEP→MPC 时 DEPDOF/INDEPDOF 直接当作 Nastran 自由度分量号
  （1–6）；若节点带 BNTRCOS 变换，MPC 应在局部系下解释，请注意核对。
- 复数载荷/复数强制位移只取实部，虚部丢弃并警告。
- 形心→剪心修正依赖 GBEAMG 的 SHCENY/SHCENZ 字段；形状卡存在但缺
  GBEAMG 时无法修正（GeniE 导出总是带 GBEAMG，手写文件才可能遇到）。
  梯形 GBARM 回退为 PBEAM 时不做剪心修正。

**约定差异（已正确处理，但值得知道）：**

- I12 = −IYZ 的符号差异源于两软件截面主轴约定不同，已按 Nastran QRG
  处理。
- 20/15/10 节点二次实体的中边节点编号按 Sesam 文档约定映射
  （IHEX：9-12 底 / 13-16 竖 / 17-20 顶），并用边端点对交叉确认。
- 实体面载荷的面号归属（SIDE 旗标 ↔ 面节点表）按 Sesam 文档插图
  整理，个别面的归属属于推断；写出前一律做外法线几何校验，
  方向不对自动翻转，但极端情况下可能加载到相邻面——重要模型请抽查。

## 构建

### 方式一：CMake（推荐，需要 MSVC / GCC / Clang）

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

产物为 `build/fem2bdf`（Windows 下 `build/Release/fem2bdf.exe`）。

### 方式二：ziglang（无需安装 Visual Studio）

适用于没有任何 C++ 工具链的 Windows 机器。zig 通过 pip 安装即可：

```bash
pip install ziglang

# Git Bash / Linux / macOS:
find extfemio/src -name '*.cpp' ! -name fem2fem.cpp ! -name bdf2bdf.cpp > srcs.txt
python -m ziglang c++ -std=c++14 -O2 -w -D_MY_CPP14_H_ -Iextfemio/include \
    @srcs.txt fem2bdf/fem2bdf.cpp -o fem2bdf.exe
```

`srcs.txt` 就是 ExtFEMIO 库的全部源文件清单（排除库自带的两个命令行
主程序）。本仓库根目录已附一份生成好的 `srcs.txt`，可直接使用。

## 用法

```bash
fem2bdf input.fem              # 生成 input.bdf
fem2bdf input.fem -o out.bdf   # 指定输出路径
```

输出为 SOL 101 完整 deck：Case Control（SUBCASE / SPC / LOAD）
+ Bulk Data。转换结束时打印统计与警告（跳过的卡片、单元类型、
缺失引用等）。

可用附带工具回读校验生成的 BDF 是否能被再次解析：

```bash
tools/bdf_check.exe out.bdf    # 输出 cards=N errors=0 即通过
```

## 图形界面（Web GUI）

`gui/` 目录提供了一个零依赖的本地 Web 界面（Python 标准库实现，
转换仍由本机 `fem2bdf.exe` 完成，文件不上传）：

```bash
python gui/server.py --open      # 默认 http://127.0.0.1:7100/
```

或直接双击 `gui/run_gui.bat`。浏览器中拖入 `.fem` 文件即可转换，
实时显示统计与警告日志，并可下载生成的 `.bdf`。

### 打包成单个 exe（免 Python 环境分发）

把 GUI、`index.html` 和 `fem2bdf.exe` 打成一个双击即用的
`fem2bdf-gui.exe`（约 10 MB），目标机器不需要装 Python：

```bash
pip install pyinstaller
python -m PyInstaller --onefile --name fem2bdf-gui --console \
    --add-data "gui/index.html;." --add-binary "fem2bdf.exe;." \
    gui/server.py
# 产物：dist/fem2bdf-gui.exe（也可用 fem2bdf-gui.spec 重新构建）
```

打包版双击启动后会自动打开浏览器（可用 `--no-open` 关闭，
`--port` 改端口）；运行时把内嵌的 `fem2bdf.exe` 解压到临时目录
调用，转换逻辑与开发版完全一致。

## 测试

```bash
python make_test_fem.py        # 重新生成 test_model.fem（167 张卡片）
./fem2bdf.exe test_model.fem -o test_model.bdf
./tools/bdf_check.exe test_model.bdf   # 回读校验，errors 应为 0
python verify_bdf.py           # 字段级逐项验证，应输出 ALL CHECKS PASSED
```

测试模型覆盖：32 个 GRID、2 个 MAT1 + 1 个 MORSMEL→MAT2、CQUAD4×2
（其中一个 GEONO 故意不定义 GELTH，触发默认厚度兜底）/ CTRIA3×2
（其中一个用 MORSMEL 各向异性材料）/ CBEAM×14（覆盖 PBEAM
与 12 种 PBEAML 截面各一条，其中两条带 BELFIX 铰接释放、一条带
GECC 局部偏心旋转；另有一条梁引用不存在的截面卡，触发跳过路径）/
CROD / CHEXA×2（含一个镜像单元自动翻转）/
CTETRA / CPENTA / CBUSH / CONM1 / CMASS2×4、2 个载荷工况、
PLOAD4×4、SPC1×2 组 + SPC×6、BLDEP→MPC 多点约束、
GSETMEMB/TDSETNAM 集合注释与 TEXT 用户文本注释。警告路径（弹簧
非对角耦合、梯形截面回退、LOTYP=2/3 跳过、随动压力、缺 BNDISPL、
GUSYI 腹板偏心、GTONP 高度约定、圆角丢弃、BELFIX 部分固定与弹性
铰、MORSMEL 的 Q 向量丢弃、缺厚度补零占位、缺截面跳过单元）
也都会被测试模型触发并出现在日志中，
方便核对警告机制本身是否工作。

## 目录结构

```
fem2bdf/
├── CMakeLists.txt        # CMake 构建脚本
├── fem2bdf/
│   └── fem2bdf.cpp       # 转换器源码（单文件）
├── extfemio/             # ExtFEMIO 库源码（vendored，含 4 处必要补丁）
├── tools/
│   └── bdf_check.cpp     # 回读校验工具（用 ExtFEMIO 重新解析生成的 BDF）
├── gui/
│   ├── server.py         # Web GUI 服务（仅标准库，调用 fem2bdf.exe）
│   ├── index.html        # GUI 页面（拖放上传 / 日志 / 下载）
│   └── run_gui.bat       # Windows 双击启动器
├── make_test_fem.py      # 生成覆盖各类卡片的测试 FEM 模型
├── verify_bdf.py         # 对生成的 BDF 做字段级逐项验证
├── test_model.fem        # 测试输入
└── test_model.bdf        # 对应转换结果（已逐项核对）
```

## 对上游库的修改（vendored 副本内，共 4 处）

1. `extfemio/include/bdf/types.h`：两个 friend operator 定义去掉 `::`
   限定（Clang/zig 不接受带限定的 friend 定义）。
2. `extfemio/src/bdf/bdf_cards_cquad4.cpp`：修复上游 bug —— 第一个
   构造函数把 `G3, G3` 传给了壳单元基类，应为 `G3, G4`。
3. 编译期定义 `_MY_CPP14_H_`（CMake 中已配置）：跳过库自带的
   `std::make_unique` C++11 回补实现，避免与 C++14 标准库冲突。
4. `extfemio/src/fem/fem_cards_gecc.cpp`：修复上游 bug —— 读取构造
   函数没有设置卡片类型标记，GECC 读出来全是 UNKNOWN，无法分发。

## 开源协议

本项目以 **BSD 0-Clause（0BSD）** 发布，见根目录 `LICENSE` 文件：
任何人可以无条件使用、复制、修改、分发本软件，商用或闭源均可，
无需保留版权声明。

内置的 ExtFEMIO 库（`extfemio/` 目录）由 DNV GL SE 以 **MIT
License** 发布，其版权声明与许可证文本保留在该目录内，
再分发时请一并保留（这是 MIT 的唯一要求）。

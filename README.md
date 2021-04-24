# PB3D

**This application is licensed under GNU GPLv3 except for the logo of the University of Canterbury (`UC.ico` and `UCBLACK.svg`).**

[![Build status](https://ci.appveyor.com/api/projects/status/udvfy5qlnfjhjqg5?svg=true)](https://ci.appveyor.com/project/TLCFEM/pb3d)

---

## How To Use

The `PB3D` is a preprocessor to prepare input file for analysis of timber structures. As in many programs, the three axes are presented by RGB colors. That is, red line represents $$x$$-axis, green line represents $$y$$-axis and blue line represents $$z$$-axis. The creation of model follows a quite straightforward flow:

1. define node or node grids,
2. define materials and sections for frame and wall elements,
3. define elements,
4. define boundary conditions and point masses,
5. configure analysis settings.

Users can use the tabs from left to right to accomplish those tasks.

### Node

Users can either create single nodes or rectangular grid of nodes under the `Node` tab. The tag/label is automatically managed so often users do not need to assign tags manually. The `Incre` column defines the spacing of grid lines in each direction while the `Repeat` column defines the number of grid lines in each direction.

The following operations are implemented to modify nodes.

#### Remove A Single Node

1. Select node tag from the list, the selected node will be highlighted.
2. Click `Apply`.

#### Remove Nodes With Fixed Interval

1. Assign the start, interval and end tags. The selected nodes will be highlighted.
    1. If start tag is left empty, only end node will be selected.
    2. If end tag is left empty, only start node will be selected.
    3. If interval is left empty, only start and end nodes will be selected.
2. Click `Apply`.

#### Remove Nodes On Given Plane/Line

1. Assign parameters $$a$$, $$b$$ and $$c$$ values. All nodes which satisfy $$ax+by+cz=0$$ will be highlighted.
2. Click `Apply`.

#### Change Coordinates

1. Assign new coordinates $$x$$, $$y$$ and $$z$$.
2. Click `Apply`.

### Section

For each section type, fill in all required parameters and click `Add` to add sections. Tags will be automatically increased but users can also assign tags manually. To remove sections, select sections from the list and click `Remove`. The section information will be printed for reference.

### Element

To add elements, users shall perform the following steps.

1. Select element type from the list, select orientation if necessary.
2. Select section from the list, which will be automatically updated depending on the element type. The section information will be printed for reference.
3. Select two nodes from the list.
4. To add multiple elements which follow a rectangular pattern, assign increment of two nodes along three axes and define how many elements should be added along each direction via `Repeat` column.

It is possible to split an element to a number of elements, to do so, select target element from the list, enter how many elements to split and click `Split`.

### BC & Load

A similar logic is used to define boundary conditions and loads/masses.

1. Select a start node from the list.
2. Define the increment of tag along each axis via `Incre X`, `Incre Y` and `Incre Z` and the number of nodes via `Repeat I`, `Repeat J` and `Repeat K`.
3. Go to `BC` or `Load/Mass` and assign desired conditions/values and click `Add`.

### Run Analysis

To perform the analysis, select the executable via `Select EXE`, save the model and click `Run`.

### Rendering Settings

It is possible to customize the rendering of model via all functionalities provided under `Settings` tab.


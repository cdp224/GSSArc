# GSSArc

Script to plot the arc of geostationary satellites as viewed from FSS earth stations in Europe.

The calculation is based on the notes found in the src directory using vector geometry. The earth is assumed to be a perfect sphere.

Take aways:

1. The observer's latitude dramatically influences the visibility and elevation of the arc.
2. The satellite's longitude does not translate into the azimuth.
3. The arc culminates on the latitude of the observer.

Tested only for latitudes of positions in Europe.

## Install

In the repository root directory:

```
python -m venv venv
venv\Scripts\activate
pip install numpy
pip install matplotlib
```

## Run

In the repository root directory:

```
venv\Scripts\activate
cd src
python GSOArcGraph.py
```

Output files are generated (e.g., pdf and png) into the ./output folder.
A window opens with the graph.

## Adjustments

### Change observers latitudes in the graph 

Modify the number of latitudes to be plotted according to your needs:

```
latitudes = [35, 40, 50, 60, 70, 80]
```

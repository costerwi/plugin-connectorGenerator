# plugin-connectorGenerator
Abaqus CAE plugin to find and create couplings and wire features connecting
circular edges matching the size and Parts of the two edges selected by the user.

This is intended to make it much more convenient to fasten layers of midsurface
together using connectors when several circular fastener holes are available in
the geometry but there are no fasteners.

1. Prompts the user to select two circular edges to connect
2. Creates reference points at center of each circular edge
3. Creates kinematic coupling from each edge to its center reference point
4. Connects center reference points using a wire connector
5. Searches for other opportunities to connect edges of the same size in the selected parts and repeates steps 2-4 above
6. Creates an assembly set for all added connector elements. These will require a connector element section assignment.

## Installation

1. Download and unzip the [latest version](https://github.com/costerwi/plugin-connectorGenerator/releases/latest)
2. Double-click the included `install.cmd` or manually copy files into your abaqus_plugins directory
3. Restart Abaqus CAE and you will find Legend King in the Visualizer plug-ins menu

## Example
The 32 connectors below were made by selecting a circular hole on the bottom instance and an adjacent circular hole in the middle instance.
Holes on the left side were skipped because there was no nearby match available.
Hole on the right side were skipped because the distance is larger than the original picked edges.
![image](https://github.com/user-attachments/assets/118fe5de-63fc-4ba1-acc2-765913edb989)

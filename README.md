# ChemPathTool

- https://ccse.lbl.gov/Research/Combustion/README_chemPathTool.htm

#Chemical Reaction Path Diagram Tool

##Introduction: Reaction Path Diagrams

The reaction path diagram tool (ChemPathTool) aids in the generation of chemical path diagrams for reacting flow systems.  More generally, ChemPathTool is used to manipulate graphically the layout of a bi-directed graph.  For reacting flow systems, the nodes of this graph represent chemical species, and the arrows (or, edges) connecting the nodes represent the flow of a conserved quantity between species, such as the transfer rate of a particular element.
In our systems, atoms are transferred among species through chemical reactions.  For example, when the reaction CH4(+M)<=>H+CH3(+M) proceeds, the impact of an arbitrary molecule (M) knocks one H atom from the CH4 molecule.  Thus, a carbon atom is "transferred" from the set of all CH4 molecules to the set of CH3.  The progress rate of this reaction (and its reverse, ie sticking an H atom onto a CH3 molecule) depends on the local state.  In a typical reacting flow calculation there are hundreds of such reactions, and the goal is to understand how the state of the fluid affects the reaction balance.  The path diagram turns out to be a good way to visualize the state of the chemical system.
   
The edges in a path diagram are obtained from a numerically simulated reacting flow using three basic steps:
A full state is obtained from the simulation.  The state might represent an instantaneous snapshot from a time-dependent simulation or a steady solution.  The state includes temperature, pressure, density and species distribution.
Forward and reverse reaction rates are computed for every relevant chemical reaction in the system using this state.
The reactions are "decomposed" into path diagram edge contributions and summed, and are thereby converted to values representing the rate that a particular atom is exchanged among chemical species.  The summing, or integration, may extend over a section of a multidimensional reacting flow domain, or may be conditioned on any parameter the investigator may find useful for the analysis.
 
The chemical species and atomic flux rates are loaded into our tool as graph nodes and directed edges, respectively.  The user then manipulates the layout of this graph into a geometry that is visually meaningful.  The resulting path diagram is written out as an Encapsulated PostScript file for printing or further manipulation, and the layout of the nodes can be saved and reused with another set of edge strengths.  By changing the input data (simulation parameters, region of integration, element of interest, etc), one can generate a series of path diagrams that are useful for identifying patterns in the chemical system.
There are a number of automated path generation tools available on the web.  But from our point of view, they all had a serious flaw that rendered them unusable.  All the tools were too smart, in that they found optimized ways to arrange the nodes based on the strength and connectivity of the edges.  Each run resulted in a different nodal layout, and any two cases were impossible to compare side by side.  (Moreover, there are standard species layouts in the combustion literature one might like to emulate).  The ChemPathTool allows (or forces, depending your point of view) one to specify the location of the nodes regardless of the edge strengths.  Since the user must now design graph layouts, we found it most simple to provide a simple customized drawing tool.  In the end, our simple concept led to a piece of software with two aims: easy creation of a nice layout, and application of that layout to multiple sets of edge strengths.
You may have picked up on one omission from the previous paragraph: how to generate the data in practice.  Like all analysis, this can be arbitrarily complicated.  In a later section, we present a simple example that will guide you in creating your own data.

## Comments and Acknowledgements

The key to the path diagram tool existence has been in its simplicity. We are not software packagers, and have only put a minimal effort into its development.  However, we strongly encourage suggestions or contributions in terms of additional features and bugs reports (email to MSDay@lbl.gov).  We hope that you find this drawing and analysis tool useful and inspiring of new analysis techniques and tools.

Also, we should mention that the path tool was developed upon ideas originally programmed by Bill Crutchfield.  Bill pointed us to Tkinter and created a first-cut implementation that has evolved into its present shape.  Also, Joe and I fleshed out our thoughts on reaction decomposition in part while talking with Dave Goodwin, the developer of the Cantera chemistry package.  Without input from these two, we may well have never built the tool.

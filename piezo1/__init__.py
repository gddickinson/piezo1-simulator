"""PIEZO1 Dynamic Structural Simulator.

An interactive, physics-driven 3D model of the PIEZO1 mechanosensitive ion
channel, built for both teaching and hypothesis generation.

Subpackages
-----------
``io``         Data acquisition and file parsing (mmCIF/PDB readers, fetchers).
``core``       Structure-of-arrays molecular model, selections, annotations.
``structure``  Superposition, cross-species numbering, hybrid model assembly,
               conformational morphing, dome geometry, pore profiling.
``physics``    Elastic network models, membrane mechanics, gating kinetics.
``analysis``   Variant mapping, contacts, pockets, conservation, docking.
``render``     OpenGL 4.1 impostor renderer (spheres, cylinders, cartoons).
``ui``         PyQt6 application shell, viewport widget and control panels.
``resources``  Curated JSON annotation data shipped with the package.
"""

__version__ = "0.1.0"
__all__ = ["__version__"]

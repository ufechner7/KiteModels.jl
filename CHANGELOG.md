## Changelog

### KiteModels v0.5.10

#### Added

-    it is now possible (and suggested) to use the DAE solver DFBDF.

This requires to add the following line to the settings.yaml file: 

    solver: "DFBDF"

The new solver is much faster (4x average, 1.8x worst case), has a lot less memory allocations (~ 50%) and is also much more stable in highly dynamic situations.

### KiteModels v0.5.8

#### Added
- new, non-allocating function `update_sys_state()`

### KiteModels v0.5.7

#### Changed
- improved performance by 10% by implementing custom norm() function for 3D vectors

### KiteModels v0.5.6

#### Fixed
- fix the method `clear!()` which failed for models with less than 6 tether segments

Simulations should work fine now for one to about 28 tether segments (no hard upper limit, but things become slow and the visualization ugly if you have too many segments).

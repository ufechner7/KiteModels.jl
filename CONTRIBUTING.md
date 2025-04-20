## License
By contributing code to KiteModels.jl, you are agreeing to release that code under the [MIT License](https://github.com/ufechner7/KiteModels.jl/blob/main/LICENSE).

## Releases
Before creating a new release, please check
- that all tests pass
- that all examples are part of the menu
- that all examples in the menu work
- update `create_sys_image2.jl` by comparing it with `create_sys_image.jl` using `meld`
- make sure that the scripts
  - test_installation
  - test_installation2 and
  - test_installation3 work
- test the installation on Linux for both Julia 1.10 and the latest stable Julia version
- test the installation on Windows after deleting the .julia folder

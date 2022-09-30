# 此项目为华为的光盘项目

## 项目构建 Makefile 

参考[链接](https://github.com/remonbonbon/makefile-example)
    
    - 环境搭建 `make install`
    - 项目编译 `make -j` (默认是在 debug 模式下，`make MODE=DEBUG -j`)
    - 通过 `compile_commands.json` 配置 `clangd` 提示 `bear make`
    - `vscode` 调试，新增 `.vscode/launch.json` 文件，需要自己先进行编译生成可执行文件 `app`
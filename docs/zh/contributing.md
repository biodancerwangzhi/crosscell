# 贡献指南

## 开发环境设置

```bash
# 克隆仓库
git clone https://github.com/biodancerwangzhi/crosscell.git
cd crosscell

# 构建 Docker 环境
docker-compose build dev

# 验证环境
docker-compose run --rm dev bash scripts/verify_environment.sh
```

## 运行测试

```bash
# 所有测试
docker-compose run --rm dev cargo test

# 特定测试
docker-compose run --rm dev cargo test test_name

# 显示输出
docker-compose run --rm dev cargo test -- --nocapture
```

## 代码风格

```bash
# 格式化代码
docker-compose run --rm dev cargo fmt

# 代码检查
docker-compose run --rm dev cargo clippy
```

## 项目结构

```
src/
├── anndata/      # AnnData (.h5ad) 读写
├── seurat/       # Seurat 对象处理
├── rds/          # RDS 格式读写
├── ir/           # 中间表示层
├── sparse/       # 稀疏矩阵操作
├── backend/      # 存储后端抽象
├── validation/   # 验证引擎
└── bin/          # CLI 入口
```

## 添加功能

1. 创建功能分支
2. 实现更改
3. 添加测试
4. 更新文档
5. 提交 Pull Request

## 测试指南

- 为新函数编写单元测试
- 为新功能添加集成测试
- 提交前确保所有测试通过

## 文档

- 更新 `docs/` 中的相关文档
- 保持英文和中文版本同步
- 使用清晰简洁的语言

## 提交信息

```
feat: 添加新功能
fix: 修复 bug
docs: 更新文档
test: 添加测试
refactor: 代码重构
```

## 有问题？

在 GitHub 上提交 issue。

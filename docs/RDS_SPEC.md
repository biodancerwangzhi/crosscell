# RDS 格式详细规范

## 文件结构

```
RDS 文件 = 文件头 + gzip(序列化数据)
```

## 1. 文件头（7 字节）

```
Byte 0-1:  魔数 "X\n" (0x58 0x0a)
Byte 2:    格式版本 (2 或 3)
Byte 3-4:  写入时的 R 版本 (major.minor)
Byte 5-6:  最小 R 版本 (major.minor)
```

## 2. 序列化数据格式

### 2.1 SEXP 类型标记（1 字节）

```rust
// 低 5 位：类型
// 高 3 位：标志位

type_byte = (flags << 5) | sexp_type

// 标志位：
const LEVELS_MASK: u8 = 0x20;    // bit 5: 有 levels（用于 factor）
const OBJECT_MASK: u8 = 0x40;    // bit 6: 是 S4 对象
const ATTR_MASK: u8 = 0x80;      // bit 7: 有属性
```

### 2.2 基础类型序列化

#### NULL (NILSXP = 0)
```
[type_byte]
```

#### Integer Vector (INTSXP = 13)
```
[type_byte]
[length: 4 bytes, little-endian]
[data: length * 4 bytes]
```

#### Real Vector (REALSXP = 14)
```
[type_byte]
[length: 4 bytes, little-endian]
[data: length * 8 bytes]
```

#### String Vector (STRSXP = 16)
```
[type_byte]
[length: 4 bytes, little-endian]
[string_1: CHARSXP]
[string_2: CHARSXP]
...
```

#### Character String (CHARSXP = 9)
```
[type_byte]
[length: 4 bytes, little-endian, 如果 length = -1 表示 NA]
[data: length bytes, UTF-8]
```

#### List (VECSXP = 19)
```
[type_byte]
[length: 4 bytes, little-endian]
[element_1: SEXP]
[element_2: SEXP]
...
```

### 2.3 属性列表

如果 ATTR_MASK 被设置，在对象数据后跟随：

```
[attributes: LISTSXP (pairlist)]
```

Pairlist 格式：
```
[type_byte: LISTSXP = 2]
[tag: SYMSXP]        // 属性名（如 "names", "class"）
[car: SEXP]          // 属性值
[cdr: LISTSXP or NILSXP]  // 下一个属性或 NULL
```

### 2.4 Symbol (SYMSXP = 1)

```
[type_byte]
[name: CHARSXP]
```

## 3. 常见属性

### names 属性
```
tag: "names"
value: STRSXP (字符串向量)
```

### class 属性
```
tag: "class"
value: STRSXP (字符串向量，如 c("data.frame"))
```

### dim 属性（矩阵/数组）
```
tag: "dim"
value: INTSXP (整数向量，如 c(nrow, ncol))
```

### row.names 属性（data.frame）
```
tag: "row.names"
value: STRSXP 或 INTSXP
```

## 4. 完整示例

### 简单整数向量 c(1, 2, 3)

```
文件头: 58 0a 03 04 03 03 05
gzip 压缩的数据:
  type_byte: 0x0d (INTSXP, 无标志)
  length: 03 00 00 00
  data: 01 00 00 00  02 00 00 00  03 00 00 00
```

### 带 names 的向量 c(a=1, b=2)

```
type_byte: 0x8d (INTSXP | ATTR_MASK)
length: 02 00 00 00
data: 01 00 00 00  02 00 00 00
attributes:
  type_byte: 0x02 (LISTSXP)
  tag: SYMSXP("names")
  car: STRSXP(c("a", "b"))
  cdr: NILSXP
```

## 5. 特殊值

### NA (Not Available)

- Integer NA: `0x80000000` (-2147483648)
- Real NA: 特殊的 NaN 位模式
- String NA: length = -1

### NaN, Inf, -Inf

- NaN: IEEE 754 NaN
- Inf: IEEE 754 正无穷
- -Inf: IEEE 754 负无穷

## 6. 引用和共享对象

RDS 格式支持对象引用以避免重复序列化：

```
REFSXP = 255
PERSISTSXP = 254
```

当遇到已序列化的对象时，使用引用而不是重新序列化。

## 7. 实现建议

### 最小可行实现（能被 R 读取）

1. 正确的文件头
2. 正确的 type_byte（包括标志位）
3. 基础类型：NILSXP, INTSXP, REALSXP, STRSXP, VECSXP
4. CHARSXP 的正确实现
5. 简单的属性支持（至少支持 names）

### 完整实现

1. 所有 SEXP 类型
2. 完整的属性系统
3. S4 对象支持
4. 引用和共享对象
5. 特殊值处理（NA, NaN, Inf）

## 8. 参考资料

- R Internals Manual: https://cran.r-project.org/doc/manuals/r-release/R-ints.html
- R 源码: src/main/serialize.c
- R 源码: src/include/Rinternals.h

## 9. 调试技巧

### 创建参考 RDS 文件

```r
# 简单向量
saveRDS(c(1L, 2L, 3L), "int_vec.rds")

# 带 names
saveRDS(c(a=1, b=2), "named_vec.rds")

# 列表
saveRDS(list(1:3, c("a", "b")), "list.rds")
```

### 查看二进制内容

```bash
# 查看文件头
hexdump -C file.rds | head -n 5

# 解压并查看
dd if=file.rds bs=1 skip=7 | gunzip | hexdump -C | head -n 20
```

### 对比输出

```bash
# 用 R 创建参考文件
Rscript -e 'saveRDS(list(1:3), "ref.rds")'

# 用 Rust 创建测试文件
cargo run -- write-rds test.rds

# 对比二进制
diff <(hexdump -C ref.rds) <(hexdump -C test.rds)
```

#!/usr/bin/env python3
"""
检查源码和测试中是否存在 mock 数据

用途：在项目发布前确保没有使用 mock 数据作为结果

检查内容：
1. Mock 库的使用（unittest.mock, pytest-mock 等）
2. 硬编码的测试数据
3. 假数据生成器（faker, factory 等）
4. 临时/示例数据
5. 跳过的测试（可能隐藏 mock）

使用方法：
    python3 scripts/check_no_mock_data.py
    python3 scripts/check_no_mock_data.py --strict  # 严格模式
    python3 scripts/check_no_mock_data.py --fix     # 自动修复（谨慎使用）
"""

import re
import sys
from pathlib import Path
from typing import List, Dict, Tuple
import argparse

# 配置
SRC_DIRS = ["src", "tests"]
EXCLUDE_DIRS = ["target", "node_modules", ".git", "__pycache__", "archive"]
EXCLUDE_FILES = ["check_no_mock_data.py"]  # 排除本脚本

# Mock 相关的模式
MOCK_PATTERNS = {
    # Python mock 库
    "mock_import": [
        r"from\s+unittest\.mock\s+import",
        r"from\s+mock\s+import",
        r"import\s+unittest\.mock",
        r"import\s+mock\b",
        r"from\s+pytest_mock\s+import",
    ],
    
    # Mock 对象使用
    "mock_usage": [
        r"@mock\.",
        r"@patch\(",
        r"Mock\(\)",
        r"MagicMock\(\)",
        r"\.mock\(",
        r"mocker\.",
    ],
    
    # 假数据生成器
    "fake_data": [
        r"from\s+faker\s+import",
        r"import\s+faker",
        r"Faker\(\)",
        r"factory\.Factory",
        r"FactoryBoy",
    ],
    
    # 硬编码测试数据（可疑模式）
    "hardcoded_data": [
        r"test_data\s*=\s*\[",
        r"mock_data\s*=",
        r"fake_data\s*=",
        r"dummy_data\s*=",
        r"sample_data\s*=\s*\{",
        r"MOCK_\w+\s*=",
        r"FAKE_\w+\s*=",
    ],
    
    # 跳过的测试（可能隐藏问题）
    "skipped_tests": [
        r"@pytest\.skip",
        r"@unittest\.skip",
        r"\.skip\(",
        r"skipIf\(",
        r"#\s*TODO.*mock",
        r"#\s*FIXME.*mock",
    ],
}

# Rust mock 模式
RUST_MOCK_PATTERNS = {
    "mock_usage": [
        r"#\[cfg\(test\)\].*mock",
        r"use\s+mockall",
        r"#\[automock\]",
        r"mock!\{",
        r"MockStruct",
    ],
}

# 允许的例外（合法使用场景）
ALLOWED_EXCEPTIONS = {
    # 测试框架本身的 mock（用于测试 mock 功能）
    "test_mock_functionality": [
        "test_mock.py",
        "test_mocking.py",
    ],
    
    # 文档和示例
    "documentation": [
        "README.md",
        "EXAMPLE.md",
        "tutorial",
    ],
}

class MockChecker:
    """Mock 数据检查器"""
    
    def __init__(self, strict_mode: bool = False):
        self.strict_mode = strict_mode
        self.issues: List[Dict] = []
        self.stats = {
            "files_scanned": 0,
            "issues_found": 0,
            "critical_issues": 0,
        }
    
    def check_file(self, filepath: Path) -> List[Dict]:
        """检查单个文件"""
        issues = []
        
        try:
            content = filepath.read_text(encoding='utf-8', errors='ignore')
            lines = content.split('\n')
            
            # 根据文件类型选择模式
            if filepath.suffix == '.rs':
                patterns = RUST_MOCK_PATTERNS
            else:
                patterns = MOCK_PATTERNS
            
            # 检查每种模式
            for category, pattern_list in patterns.items():
                for pattern in pattern_list:
                    for line_num, line in enumerate(lines, 1):
                        if re.search(pattern, line, re.IGNORECASE):
                            # 检查是否是允许的例外
                            if self._is_allowed_exception(filepath, line):
                                continue
                            
                            issues.append({
                                "file": str(filepath),
                                "line": line_num,
                                "category": category,
                                "pattern": pattern,
                                "content": line.strip(),
                                "severity": self._get_severity(category),
                            })
            
        except Exception as e:
            print(f"⚠️  无法读取文件 {filepath}: {e}")
        
        return issues
    
    def _is_allowed_exception(self, filepath: Path, line: str) -> bool:
        """检查是否是允许的例外"""
        # 检查文件名例外
        for exception_file in ALLOWED_EXCEPTIONS["test_mock_functionality"]:
            if exception_file in str(filepath):
                return True
        
        # 检查注释（如果是注释中的 mock，可能只是说明）
        if line.strip().startswith('#') or line.strip().startswith('//'):
            # 但如果是 TODO/FIXME 注释，仍然标记
            if 'TODO' in line or 'FIXME' in line:
                return False
            return True
        
        return False
    
    def _get_severity(self, category: str) -> str:
        """获取问题严重程度"""
        critical_categories = ["mock_import", "mock_usage", "fake_data"]
        
        if category in critical_categories:
            return "CRITICAL"
        elif category == "hardcoded_data":
            return "WARNING"
        else:
            return "INFO"
    
    def scan_directory(self, directory: Path):
        """扫描目录"""
        print(f"\n🔍 扫描目录: {directory}")
        
        for filepath in directory.rglob("*"):
            # 跳过目录
            if filepath.is_dir():
                continue
            
            # 跳过排除的目录
            if any(excluded in filepath.parts for excluded in EXCLUDE_DIRS):
                continue
            
            # 跳过排除的文件
            if filepath.name in EXCLUDE_FILES:
                continue
            
            # 只检查源码文件
            if filepath.suffix not in ['.py', '.rs', '.R']:
                continue
            
            self.stats["files_scanned"] += 1
            
            # 检查文件
            file_issues = self.check_file(filepath)
            if file_issues:
                self.issues.extend(file_issues)
                self.stats["issues_found"] += len(file_issues)
                
                # 统计严重问题
                critical = sum(1 for issue in file_issues if issue["severity"] == "CRITICAL")
                self.stats["critical_issues"] += critical
    
    def print_report(self):
        """打印检查报告"""
        print("\n" + "="*70)
        print("📊 Mock 数据检查报告")
        print("="*70)
        
        print(f"\n📁 扫描统计:")
        print(f"  - 扫描文件数: {self.stats['files_scanned']}")
        print(f"  - 发现问题数: {self.stats['issues_found']}")
        print(f"  - 严重问题数: {self.stats['critical_issues']}")
        
        if not self.issues:
            print("\n✅ 太好了！没有发现 mock 数据使用")
            return True
        
        # 按严重程度分组
        critical_issues = [i for i in self.issues if i["severity"] == "CRITICAL"]
        warning_issues = [i for i in self.issues if i["severity"] == "WARNING"]
        info_issues = [i for i in self.issues if i["severity"] == "INFO"]
        
        # 打印严重问题
        if critical_issues:
            print("\n" + "="*70)
            print("🚨 严重问题 (CRITICAL) - 必须修复")
            print("="*70)
            self._print_issues(critical_issues)
        
        # 打印警告
        if warning_issues:
            print("\n" + "="*70)
            print("⚠️  警告 (WARNING) - 建议检查")
            print("="*70)
            self._print_issues(warning_issues)
        
        # 打印信息（仅在严格模式）
        if info_issues and self.strict_mode:
            print("\n" + "="*70)
            print("ℹ️  信息 (INFO) - 可能需要注意")
            print("="*70)
            self._print_issues(info_issues)
        
        # 总结
        print("\n" + "="*70)
        print("📋 修复建议")
        print("="*70)
        print("\n1. 移除所有 mock 库的导入和使用")
        print("2. 使用真实数据进行测试")
        print("3. 如果必须使用测试数据，从文件加载而不是硬编码")
        print("4. 取消跳过的测试，确保所有测试都能通过")
        print("5. 移除所有假数据生成器")
        
        return False
    
    def _print_issues(self, issues: List[Dict]):
        """打印问题列表"""
        # 按文件分组
        by_file = {}
        for issue in issues:
            file = issue["file"]
            if file not in by_file:
                by_file[file] = []
            by_file[file].append(issue)
        
        # 打印每个文件的问题
        for file, file_issues in sorted(by_file.items()):
            print(f"\n📄 {file}")
            for issue in file_issues:
                print(f"  ├─ 行 {issue['line']}: [{issue['category']}]")
                print(f"  │  {issue['content'][:80]}")
                if len(issue['content']) > 80:
                    print(f"  │  ...")
    
    def generate_fix_suggestions(self) -> List[str]:
        """生成修复建议"""
        suggestions = []
        
        # 按文件分组问题
        by_file = {}
        for issue in self.issues:
            if issue["severity"] == "CRITICAL":
                file = issue["file"]
                if file not in by_file:
                    by_file[file] = []
                by_file[file].append(issue)
        
        for file, issues in by_file.items():
            suggestions.append(f"\n# 修复 {file}")
            suggestions.append(f"# 发现 {len(issues)} 个严重问题")
            
            # 根据问题类型给出建议
            for issue in issues:
                if "mock_import" in issue["category"]:
                    suggestions.append(f"# 行 {issue['line']}: 移除 mock 导入")
                    suggestions.append(f"# 建议: 使用真实数据或从文件加载测试数据")
                
                elif "mock_usage" in issue["category"]:
                    suggestions.append(f"# 行 {issue['line']}: 移除 mock 对象使用")
                    suggestions.append(f"# 建议: 重构测试以使用真实对象")
                
                elif "fake_data" in issue["category"]:
                    suggestions.append(f"# 行 {issue['line']}: 移除假数据生成器")
                    suggestions.append(f"# 建议: 使用真实数据集")
        
        return suggestions

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="检查源码中是否存在 mock 数据",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python3 scripts/check_no_mock_data.py              # 标准检查
  python3 scripts/check_no_mock_data.py --strict     # 严格模式
  python3 scripts/check_no_mock_data.py --fix        # 生成修复建议
        """
    )
    
    parser.add_argument(
        "--strict",
        action="store_true",
        help="严格模式：报告所有可疑模式"
    )
    
    parser.add_argument(
        "--fix",
        action="store_true",
        help="生成修复建议"
    )
    
    parser.add_argument(
        "--output",
        type=str,
        help="将报告输出到文件"
    )
    
    args = parser.parse_args()
    
    # 创建检查器
    checker = MockChecker(strict_mode=args.strict)
    
    # 扫描源码目录
    project_root = Path.cwd()
    for src_dir in SRC_DIRS:
        dir_path = project_root / src_dir
        if dir_path.exists():
            checker.scan_directory(dir_path)
    
    # 打印报告
    success = checker.print_report()
    
    # 生成修复建议
    if args.fix and checker.issues:
        print("\n" + "="*70)
        print("🔧 修复建议")
        print("="*70)
        suggestions = checker.generate_fix_suggestions()
        for suggestion in suggestions:
            print(suggestion)
    
    # 输出到文件
    if args.output:
        output_path = Path(args.output)
        with output_path.open('w', encoding='utf-8') as f:
            f.write("Mock 数据检查报告\n")
            f.write("="*70 + "\n\n")
            f.write(f"扫描文件数: {checker.stats['files_scanned']}\n")
            f.write(f"发现问题数: {checker.stats['issues_found']}\n")
            f.write(f"严重问题数: {checker.stats['critical_issues']}\n\n")
            
            for issue in checker.issues:
                f.write(f"{issue['file']}:{issue['line']} [{issue['severity']}]\n")
                f.write(f"  {issue['content']}\n\n")
        
        print(f"\n📝 报告已保存到: {output_path}")
    
    # 返回退出码
    if checker.stats["critical_issues"] > 0:
        print("\n❌ 发现严重问题，请修复后再发布")
        return 1
    elif checker.stats["issues_found"] > 0 and args.strict:
        print("\n⚠️  发现警告，建议检查")
        return 1
    else:
        print("\n✅ 检查通过，可以发布")
        return 0

if __name__ == "__main__":
    sys.exit(main())

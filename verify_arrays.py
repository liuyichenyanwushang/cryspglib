#!/usr/bin/env python3
"""
验证C和Rust代码中的数组一致性
"""

import re

def extract_c_array(filename, array_name):
    """从C文件中提取数组定义"""
    with open(filename, 'r') as f:
        content = f.read()

    # 查找数组定义
    pattern = rf'static\s+int\s+{array_name}\[(\d+)\]\s*=\s*\{{(.*?)\}};'
    match = re.search(pattern, content, re.DOTALL)
    if not match:
        return None

    size = int(match.group(1))
    body = match.group(2)

    # 提取所有数字
    numbers = re.findall(r'(-?\d+)', body)
    values = [int(n) for n in numbers]

    return size, values

def extract_rust_array(filename, array_name):
    """从Rust文件中提取数组定义"""
    with open(filename, 'r') as f:
        content = f.read()

    # 查找数组定义
    pattern = rf'static\s+{array_name}:\s*\[i32;\s*(\d+)\]\s*=\s*\[(.*?)\];'
    match = re.search(pattern, content, re.DOTALL)
    if not match:
        return None

    size = int(match.group(1))
    body = match.group(2)

    # 提取所有数字
    numbers = re.findall(r'(-?\d+)', body)
    values = [int(n) for n in numbers]

    return size, values

def main():
    print("验证数组一致性")
    print("=" * 60)

    # 检查 SPACEGROUP_TO_HALL_NUMBER
    print("\n1. SPACEGROUP_TO_HALL_NUMBER:")
    c_result = extract_c_array('src_c/spacegroup.c', 'spacegroup_to_hall_number')
    rust_result = extract_rust_array('src/spacegroup.rs', 'SPACEGROUP_TO_HALL_NUMBER')

    if c_result and rust_result:
        c_size, c_values = c_result
        r_size, r_values = rust_result

        print(f"   C: size={c_size}, values={len(c_values)}")
        print(f"   Rust: size={r_size}, values={len(r_values)}")

        if c_size == r_size and len(c_values) == len(r_values) and c_values == r_values:
            print("   ✅ 数组一致")
        else:
            print("   ❌ 数组不一致")
            if c_size != r_size:
                print(f"     大小不匹配: C={c_size}, Rust={r_size}")
            if len(c_values) != len(r_values):
                print(f"     值数量不匹配: C={len(c_values)}, Rust={len(r_values)}")
            if c_values != r_values:
                print(f"     值不匹配")
                for i, (c_val, r_val) in enumerate(zip(c_values, r_values)):
                    if c_val != r_val:
                        print(f"     索引 {i}: C={c_val}, Rust={r_val}")
    else:
        print("   ⚠ 无法提取数组")

    # 检查 LAYER_GROUP_TO_HALL_NUMBER
    print("\n2. LAYER_GROUP_TO_HALL_NUMBER:")
    c_result = extract_c_array('src_c/spacegroup.c', 'layer_group_to_hall_number')
    rust_result = extract_rust_array('src/spacegroup.rs', 'LAYER_GROUP_TO_HALL_NUMBER')

    if c_result and rust_result:
        c_size, c_values = c_result
        r_size, r_values = rust_result

        print(f"   C: size={c_size}, values={len(c_values)}")
        print(f"   Rust: size={r_size}, values={len(r_values)}")

        if c_size == r_size and len(c_values) == len(r_values) and c_values == r_values:
            print("   ✅ 数组一致")
        else:
            print("   ❌ 数组不一致")
            if c_size != r_size:
                print(f"     大小不匹配: C={c_size}, Rust={r_size}")
            if len(c_values) != len(r_values):
                print(f"     值数量不匹配: C={len(c_values)}, Rust={len(r_values)}")

            # 检查是否C数组声明大小大于实际值
            if c_size > len(c_values):
                print(f"     ⚠ C数组声明大小{c_size}大于实际值数量{len(c_values)}")
                print(f"     这可能是一个bug：数组未完全初始化")

    print("\n" + "=" * 60)

if __name__ == '__main__':
    main()
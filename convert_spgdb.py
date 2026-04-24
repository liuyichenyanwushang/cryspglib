#!/usr/bin/env python3
"""
Convert spg_database.c static data to Rust static arrays.
Usage: ./convert_spgdb.py spg_database.c spg_database.rs
"""

import re
import sys

CENTERING_MAP = {
    'CENTERING_ERROR': 'Centering::Error',
    'PRIMITIVE': 'Centering::Primitive',
    'BODY': 'Centering::Body',
    'FACE': 'Centering::Face',
    'A_FACE': 'Centering::AFace',
    'B_FACE': 'Centering::BFace',
    'C_FACE': 'Centering::CFace',
    'R_CENTER': 'Centering::RCenter',
}

def remove_c_comments(text):
    """移除 C 风格的注释 /* ... */"""
    return re.sub(r'/\*.*?\*/', '', text, flags=re.DOTALL)

def extract_array_body(content, array_name):
    """
    通过查找数组名和大括号直接提取数组内容，避免正则表达式的复杂匹配。
    """
    idx = content.find(array_name)
    if idx == -1:
        raise ValueError(f"未找到数组 {array_name}")
    # 找到等号
    eq_idx = content.find('=', idx)
    if eq_idx == -1:
        raise ValueError(f"数组 {array_name} 后未找到等号")
    # 找到第一个 {
    brace_start = content.find('{', eq_idx)
    if brace_start == -1:
        raise ValueError(f"数组 {array_name} 后未找到左大括号")
    # 匹配对应的 }
    brace_level = 0
    brace_end = -1
    for i in range(brace_start, len(content)):
        if content[i] == '{':
            brace_level += 1
        elif content[i] == '}':
            brace_level -= 1
            if brace_level == 0:
                brace_end = i
                break
    if brace_end == -1:
        raise ValueError(f"数组 {array_name} 的右大括号未找到")
    # 提取内容，不包括外层大括号
    return content[brace_start+1:brace_end].strip()

def split_struct_initializers(struct_array_body):
    """将结构体数组内容拆分为单个结构体初始化字符串"""
    structs = []
    brace_level = 0
    start = 0
    for i, ch in enumerate(struct_array_body):
        if ch == '{':
            if brace_level == 0:
                start = i + 1
            brace_level += 1
        elif ch == '}':
            brace_level -= 1
            if brace_level == 0:
                structs.append(struct_array_body[start:i].strip())
    return structs

def parse_struct_fields(struct_body):
    """解析结构体内部的字段"""
    fields = []
    current = ''
    in_string = False
    for ch in struct_body:
        if ch == '"':
            in_string = not in_string
            current += ch
        elif ch == ',' and not in_string:
            fields.append(current.strip())
            current = ''
        else:
            current += ch
    if current:
        fields.append(current.strip())
    return fields

def clean_string(s):
    """去除 C 字符串字面量的引号"""
    if s.startswith('"') and s.endswith('"'):
        return s[1:-1]
    return s

def parse_spacegroup_type(struct_body):
    """解析 SpacegroupType 结构体，生成 Rust 初始化代码"""
    fields = parse_struct_fields(struct_body)
    if len(fields) != 9:
        raise ValueError(f"结构体字段数量错误: 预期 9，实际 {len(fields)}: {fields}")

    number = int(fields[0])
    schoenflies = clean_string(fields[1])
    hall_symbol = clean_string(fields[2])
    international = clean_string(fields[3])
    international_full = clean_string(fields[4])
    international_short = clean_string(fields[5])
    choice = clean_string(fields[6])
    centering_c = fields[7].strip()
    pointgroup_number = int(fields[8])

    centering = CENTERING_MAP.get(centering_c, f"/* 未知 {centering_c} */")
    return (f"    RawSpacegroupType {{\n"
            f"        number: {number},\n"
            f"        schoenflies: \"{schoenflies}\",\n"
            f"        hall_symbol: \"{hall_symbol}\",\n"
            f"        international: \"{international}\",\n"
            f"        international_full: \"{international_full}\",\n"
            f"        international_short: \"{international_short}\",\n"
            f"        choice: \"{choice}\",\n"
            f"        centering: {centering},\n"
            f"        pointgroup_number: {pointgroup_number},\n"
            f"    }}")

def parse_int_array(body):
    """从逗号分隔的整数列表中提取所有整数"""
    return [int(x) for x in re.findall(r'-?\d+', body)]

def parse_pair_array(body):
    """解析二维数组中的每一对 {x, y}"""
    pairs = []
    for m in re.finditer(r'{\s*(\d+)\s*,\s*(\d+)\s*}', body):
        pairs.append((int(m.group(1)), int(m.group(2))))
    return pairs

def generate_rust(c_code):
    rust_lines = []

    # 移除注释，避免干扰解析
    c_code = remove_c_comments(c_code)

    rust_lines.append("//! Automatically generated from spg_database.c")
    rust_lines.append("//! DO NOT EDIT MANUALLY\n")
    rust_lines.append("use crate::mathfunc::{Mat3I, Vec3};\n")
    rust_lines.append("#[derive(Debug, Clone, Copy, PartialEq)]")
    rust_lines.append("pub enum Centering {")
    rust_lines.append("    Error,")
    rust_lines.append("    Primitive,")
    rust_lines.append("    Body,")
    rust_lines.append("    Face,")
    rust_lines.append("    AFace,")
    rust_lines.append("    BFace,")
    rust_lines.append("    CFace,")
    rust_lines.append("    RCenter,")
    rust_lines.append("}\n")
    rust_lines.append("/// Raw space group type data as stored in the database.")
    rust_lines.append("#[derive(Debug, Clone, Copy)]")
    rust_lines.append("pub struct RawSpacegroupType {")
    rust_lines.append("    pub number: i32,")
    rust_lines.append("    pub schoenflies: &'static str,")
    rust_lines.append("    pub hall_symbol: &'static str,")
    rust_lines.append("    pub international: &'static str,")
    rust_lines.append("    pub international_full: &'static str,")
    rust_lines.append("    pub international_short: &'static str,")
    rust_lines.append("    pub choice: &'static str,")
    rust_lines.append("    pub centering: Centering,")
    rust_lines.append("    pub pointgroup_number: i32,")
    rust_lines.append("}\n")

    # ---------- spacegroup_types ----------
    body = extract_array_body(c_code, 'spacegroup_types')
    structs = split_struct_initializers(body)
    rust_lines.append(f"pub static SPACEGROUP_TYPES: [RawSpacegroupType; {len(structs)}] = [")
    for s in structs:
        rust_lines.append(parse_spacegroup_type(s))
    rust_lines.append("];\n")

    # ---------- layer_group_types ----------
    body = extract_array_body(c_code, 'layer_group_types')
    structs = split_struct_initializers(body)
    rust_lines.append(f"pub static LAYER_GROUP_TYPES: [RawSpacegroupType; {len(structs)}] = [")
    for s in structs:
        rust_lines.append(parse_spacegroup_type(s))
    rust_lines.append("];\n")

    # ---------- symmetry_operations ----------
    body = extract_array_body(c_code, 'symmetry_operations')
    ints = parse_int_array(body)
    rust_lines.append(f"pub static SYMMETRY_OPERATIONS: [i32; {len(ints)}] = [")
    for i in range(0, len(ints), 8):
        chunk = ints[i:i+8]
        rust_lines.append("    " + ", ".join(str(x) for x in chunk) + ",")
    rust_lines.append("];\n")

    # ---------- symmetry_operation_index ----------
    body = extract_array_body(c_code, 'symmetry_operation_index')
    pairs = parse_pair_array(body)
    rust_lines.append(f"pub static SYMMETRY_OPERATION_INDEX: [[i32; 2]; {len(pairs)}] = [")
    for a, b in pairs:
        rust_lines.append(f"    [{a}, {b}],")
    rust_lines.append("];\n")

    # ---------- layer_symmetry_operation_index ----------
    body = extract_array_body(c_code, 'layer_symmetry_operation_index')
    pairs = parse_pair_array(body)
    rust_lines.append(f"pub static LAYER_SYMMETRY_OPERATION_INDEX: [[i32; 2]; {len(pairs)}] = [")
    for a, b in pairs:
        rust_lines.append(f"    [{a}, {b}],")
    rust_lines.append("];\n")

    return "\n".join(rust_lines)

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.c output.rs")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    with open(input_file, 'r', encoding='utf-8') as f:
        c_code = f.read()

    rust_code = generate_rust(c_code)

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(rust_code)

    print(f"成功生成 {output_file}")

if __name__ == "__main__":
    main()

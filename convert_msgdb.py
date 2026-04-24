#!/usr/bin/env python3
"""
Convert msg_database.c static data to Rust static arrays.

Usage: ./convert_msgdb.py src_c/msg_database.txt src/msg_database.rs
"""

import re
import sys

def remove_c_comments(text):
    return re.sub(r'/\*.*?\*/', '', text, flags=re.DOTALL)

def extract_array_body(content, array_name):
    """通过查找数组名和大括号直接提取数组内容"""
    idx = content.find(array_name)
    if idx == -1:
        raise ValueError(f"未找到数组 {array_name}")
    eq_idx = content.find('=', idx)
    if eq_idx == -1:
        raise ValueError(f"数组 {array_name} 后未找到等号")
    brace_start = content.find('{', eq_idx)
    if brace_start == -1:
        raise ValueError("未找到左大括号")
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
        raise ValueError("右大括号未找到")
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
    """解析结构体内部的字段（支持字符串中的逗号）"""
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
    if s.startswith('"') and s.endswith('"'):
        return s[1:-1]
    return s

def parse_magnetic_spacegroup_type(struct_body):
    """解析 MagneticSpacegroupType 结构体"""
    fields = parse_struct_fields(struct_body)
    if len(fields) != 6:
        raise ValueError(f"字段数量错误: 预期 6，实际 {len(fields)}: {fields}")
    uni_number = int(fields[0])
    litvin_number = int(fields[1])
    bns_number = clean_string(fields[2])
    og_number = clean_string(fields[3])
    number = int(fields[4])
    mtype = int(fields[5])
    return (f"    MagneticSpacegroupType {{\n"
            f"        uni_number: {uni_number},\n"
            f"        litvin_number: {litvin_number},\n"
            f"        bns_number: \"{bns_number}\",\n"
            f"        og_number: \"{og_number}\",\n"
            f"        number: {number},\n"
            f"        type_: {mtype},\n"
            f"    }},")

def parse_int_array(body):
    """提取所有整数"""
    return [int(x) for x in re.findall(r'-?\d+', body)]

def parse_pair_array(body):
    """解析 {x, y} 对"""
    pairs = []
    for m in re.finditer(r'{\s*(-?\d+)\s*,\s*(-?\d+)\s*}', body):
        pairs.append((int(m.group(1)), int(m.group(2))))
    return pairs

def split_nested_initializers(body_or_file, inner_count):
    """Split an array body into top-level entries using regex.
    Each entry is '{' pairs '}' where pairs are '{a,b}, {c,d}, ...'.
    Returns list of list of (a,b) tuples."""
    # Use regex to find top-level { ... } entries correctly
    # Strategy: find balanced braces at level 1
    brace_level = 0
    entries = []
    start = None
    for i, ch in enumerate(body_or_file):
        if ch == '{':
            if brace_level == 0:
                start = i + 1
            brace_level += 1
        elif ch == '}':
            brace_level -= 1
            if brace_level == 0 and start is not None:
                entry_content = body_or_file[start:i].strip()
                # Now find all {a,b} pairs within entry_content
                pairs = []
                pb = 0
                ps = None
                for j, c2 in enumerate(entry_content):
                    if c2 == '{':
                        if pb == 0:
                            ps = j + 1
                        pb += 1
                    elif c2 == '}':
                        pb -= 1
                        if pb == 0 and ps is not None:
                            text = entry_content[ps:j].strip()
                            vals = [int(x) for x in text.split(',')]
                            if len(vals) >= 2:
                                pairs.append((vals[0], vals[1]))
                            ps = None
                entries.append(pairs)
                start = None
    return entries

def split_nested_seven(body):
    """Split [][18][7] array into entries, each with up to 18 {7-int} groups."""
    entries = []
    brace_level = 0
    entry_start = None
    for i, ch in enumerate(body):
        if ch == '{':
            brace_level += 1
            if brace_level == 1:
                entry_start = i + 1
        elif ch == '}':
            brace_level -= 1
            if brace_level == 0 and entry_start is not None:
                entry_body = body[entry_start:i].strip()
                # Parse {a,b,c,d,e,f,g} groups from entry body
                groups = []
                gb_level = 0
                gb_start = None
                for j, c2 in enumerate(entry_body):
                    if c2 == '{':
                        if gb_level == 0:
                            gb_start = j + 1
                        gb_level += 1
                    elif c2 == '}':
                        gb_level -= 1
                        if gb_level == 0 and gb_start is not None:
                            group_text = entry_body[gb_start:j].strip()
                            vals = [int(x) for x in group_text.split(',') if x.strip()]
                            if len(vals) >= 7:
                                groups.append(vals[:7])
                            gb_start = None
                entries.append(groups)
                start = None
    return entries

def generate_rust(c_code):
    rust_lines = []
    c_code = remove_c_comments(c_code)

    rust_lines.append("// Automatically generated from msg_database.txt")
    rust_lines.append("// DO NOT EDIT MANUALLY\n")

    # ---------- MagneticSpacegroupType struct ----------
    rust_lines.append("#[derive(Debug, Clone, Copy)]")
    rust_lines.append("pub struct MagneticSpacegroupType {")
    rust_lines.append("    pub uni_number: i32,")
    rust_lines.append("    pub litvin_number: i32,")
    rust_lines.append("    pub bns_number: &'static str,")
    rust_lines.append("    pub og_number: &'static str,")
    rust_lines.append("    pub number: i32,")
    rust_lines.append("    pub type_: i32,")  # 'type' is keyword
    rust_lines.append("}\n")

    # ---------- magnetic_spacegroup_types ----------
    body = extract_array_body(c_code, 'magnetic_spacegroup_types')
    structs = split_struct_initializers(body)
    rust_lines.append(f"pub static MAGNETIC_SPACEGROUP_TYPES: [MagneticSpacegroupType; {len(structs)}] = [")
    for s in structs:
        rust_lines.append(parse_magnetic_spacegroup_type(s))
    rust_lines.append("];\n")

    # ---------- magnetic_spacegroup_hall_mapping ----------
    body = extract_array_body(c_code, 'magnetic_spacegroup_hall_mapping')
    pairs = parse_pair_array(body)
    rust_lines.append(f"pub static MAGNETIC_SPACEGROUP_HALL_MAPPING: [[i32; 2]; {len(pairs)}] = [")
    for a, b in pairs:
        rust_lines.append(f"    [{a}, {b}],")
    rust_lines.append("];\n")

    # ---------- magnetic_spacegroup_uni_mapping ----------
    body = extract_array_body(c_code, 'magnetic_spacegroup_uni_mapping')
    pairs = parse_pair_array(body)
    rust_lines.append(f"pub static MAGNETIC_SPACEGROUP_UNI_MAPPING: [[i32; 2]; {len(pairs)}] = [")
    for a, b in pairs:
        rust_lines.append(f"    [{a}, {b}],")
    rust_lines.append("];\n")

    # ---------- magnetic_spacegroup_operation_index ----------
    body = extract_array_body(c_code, 'magnetic_spacegroup_operation_index')
    # Format: [][18][2] - each top-level group has variable pairs, pad to 18
    entries = split_nested_initializers(body, 2)
    rust_lines.append(f"pub static MAGNETIC_SPACEGROUP_OPERATION_INDEX: [[[i32; 2]; 18]; {len(entries)}] = [")
    for pairs in entries:
        rows = []
        for p in pairs:
            rows.append(f"[{p[0]}, {p[1]}]")
        # Pad to 18 pairs
        while len(rows) < 18:
            rows.append("[0, 0]")
        indent = "        "
        rust_lines.append("[" + ", ".join(rows) + "],")
    rust_lines.append("];\n")

    # ---------- magnetic_symmetry_operations ----------
    body = extract_array_body(c_code, 'magnetic_symmetry_operations')
    ints = parse_int_array(body)
    rust_lines.append(f"pub static MAGNETIC_SYMMETRY_OPERATIONS: [i32; {len(ints)}] = [")
    for i in range(0, len(ints), 8):
        chunk = ints[i:i+8]
        rust_lines.append("    " + ", ".join(str(x) for x in chunk) + ",")
    rust_lines.append("];\n")

    # ---------- alternative_transformations ----------
    body = extract_array_body(c_code, 'alternative_transformations')
    # Format: [][18][7] - each top-level group has variable {a,b,c,d,e,f,g} sub-groups (7 ints)
    alt_entries = split_nested_seven(body)
    rust_lines.append(f"pub static ALTERNATIVE_TRANSFORMATIONS: [[[i32; 7]; 18]; {len(alt_entries)}] = [")
    for groups in alt_entries:
        rows = []
        for g in groups:
            rows.append("[" + ", ".join(str(x) for x in g) + "]")
        while len(rows) < 18:
            rows.append("[0, 0, 0, 0, 0, 0, 0]")
        indent = "        "
        rust_lines.append("[" + ", ".join(rows) + "],")
    rust_lines.append("];\n")

    # ---------- msgdb_get_uni_candidates ----------
    # Already covered by MAGNETIC_SPACEGROUP_HALL_MAPPING
    # We'll put tiny logic functions in msg_database.rs manually

    return "\n".join(rust_lines)

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.c output.rs")
        sys.exit(1)

    with open(sys.argv[1], 'r', encoding='utf-8') as f:
        c_code = f.read()

    rust_code = generate_rust(c_code)

    with open(sys.argv[2], 'w', encoding='utf-8') as f:
        f.write(rust_code)

    print(f"成功生成 {sys.argv[2]}")

if __name__ == "__main__":
    main()

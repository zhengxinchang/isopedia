#!/usr/bin/env python3

import random
from tqdm import tqdm
from collections import defaultdict


def calculate_increment(prev_offset: int, prev_length: int) -> int:
    """
    基于前一条记录的offset和length计算增量
    
    参数:
        prev_offset: 前一条记录的offset值
        prev_length: 前一条记录的length值
    
    返回:
        增量值，确保新的length在50-400范围内
    """
    # 使用前一条记录的信息计算基础增量
    base_increment = (prev_offset + prev_length) % 20 - 10  # 产生-10到+10的增量
    return base_increment


def generate_mock_data(filename: str, num_records: int = 100000, duplicate_ratio: float = 0.3):
    """
    生成模拟的基因组数据，包含重复的坐标和基于前一条记录计算的长度

    参数:
        filename: 输出文件名
        num_records: 记录总数
        duplicate_ratio: 重复坐标的比例
    """
    with open(filename, 'w') as f:
        records = []
        # 生成坐标字典，用于追踪重复
        coord_count = defaultdict(int)

        # 初始化第一条记录的length
        prev_length = random.randint(50, 400)
        prev_offset = 0

        # 先生成基础记录
        base_records = []
        num_unique = int(num_records * (1 - duplicate_ratio))

        for offset in tqdm(range(num_unique), desc="Generating unique records"):
            chr_num = 0
            pos = random.randint(0, num_records)
            
            # 计算新的length
            increment = calculate_increment(prev_offset, prev_length)
            new_length = prev_length + increment
            
            # 确保length在有效范围内
            new_length = max(50, min(400, new_length))
            
            coord = (chr_num, pos)
            coord_count[coord] = 1
            base_records.append((chr_num, pos, offset, new_length))
            
            # 更新前一条记录的值
            prev_length = new_length
            prev_offset = offset

        # 生成重复记录
        for offset in tqdm(range(num_unique, num_records), desc="Generating duplicate records"):
            # 随机选择一个已存在的坐标进行重复
            if base_records:  # 确保有基础记录可选
                chr_num, pos = random.choice(base_records)[:2]
                
                # 计算新的length
                increment = calculate_increment(prev_offset, prev_length)
                new_length = prev_length + increment
                
                # 确保length在有效范围内
                new_length = max(50, min(400, new_length))
                
                coord_count[(chr_num, pos)] += 1
                records.append((chr_num, pos, offset, new_length))
                
                # 更新前一条记录的值
                prev_length = new_length
                prev_offset = offset

        # 合并基础记录和重复记录
        records = base_records + records

        # 按染色体号和位置排序
        records.sort(key=lambda x: (x[0], x[1]))

        # 写入文件
        for chr_num, pos, offset, length in records:
            f.write(f"{chr_num}:{pos}:{offset}:{length}\n")

        # 打印重复统计信息
        print("\n重复坐标统计:")
        duplicate_coords = [(coord, count)
                          for coord, count in coord_count.items() if count > 1]
        print(f"重复坐标数量: {len(duplicate_coords)}")
        print(f"最高重复次数: {max(coord_count.values()) if coord_count else 0}")


def verify_data(filename: str):
    """验证生成的数据，包括重复检查和length的连续性验证"""
    with open(filename, 'r') as f:
        lines = f.readlines()
        print(f"总行数: {len(lines)}")

        # 检查几行数据
        print("\n示例数据:")
        for line in lines[:5]:
            print(line.strip())

        # 验证length的连续性
        print("\n验证length的连续性:")
        prev_offset = None
        prev_length = None
        for i, line in enumerate(lines[:10]):  # 显示前10条记录的验证
            parts = line.strip().split(':')
            curr_offset = int(parts[2])
            curr_length = int(parts[3])
            
            if prev_offset is not None:
                expected_increment = calculate_increment(prev_offset, prev_length)
                expected_length = max(50, min(400, prev_length + expected_increment))
                print(f"记录 {i}: offset={curr_offset}, length={curr_length}, " 
                      f"基于前值计算增量={expected_increment}")
            
            prev_offset = curr_offset
            prev_length = curr_length

        # 统计重复坐标
        coord_dict = defaultdict(list)
        for line in lines:
            parts = line.strip().split(':')
            chr_num, pos, offset, length = map(int, parts)
            coord_dict[f"{chr_num}:{pos}"].append((offset, length))

        # 显示重复情况
        duplicates = {k: v for k, v in coord_dict.items() if len(v) > 1}
        print("\n重复坐标示例:")
        for k, v in list(duplicates.items())[:5]:
            print(f"坐标 {k} 出现 {len(v)} 次，记录: {v}")


if __name__ == "__main__":
    output_file = "./test/mock_interim.txt"
    num_records = 100000000
    duplicate_ratio = 0.3  # 30%的坐标会重复
    generate_mock_data(output_file, num_records, duplicate_ratio)
    verify_data(output_file)
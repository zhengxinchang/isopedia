import sqlite3
import random
import time
from typing import List, Tuple

def perform_random_queries(db_file: str = "test/sqlite3.db", 
                         num_queries: int = 10000000) -> Tuple[float, List[Tuple[float, int, List[Tuple[int, ...]]]]]:
    """
    执行随机查询测试
    返回：总用时和详细查询记录
    """
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    
    # 获取col2的最大值和最小值
    cursor.execute('SELECT MIN(col2), MAX(col2) FROM data_table')
    min_val, max_val = cursor.fetchone()
    
    query_records = []  # 存储每次查询的记录
    total_start_time = time.time()
    
    try:
        for i in range(num_queries):
            # 生成随机值
            random_val = random.randint(min_val, max_val)
            
            # 执行查询并记录时间
            query_start_time = time.time()
            cursor.execute('SELECT * FROM data_table WHERE col2 = ?', (random_val,))
            results = cursor.fetchall()
            query_time = time.time() - query_start_time
            
            # 存储查询记录
            query_records.append((query_time, random_val, results))
            
            # 每100000次查询显示进度
            if (i + 1) % 100000 == 0:
                print(f"已完成 {i + 1} 次查询")
        
        total_time = time.time() - total_start_time
        return total_time, query_records
        
    finally:
        cursor.close()
        conn.close()

def analyze_results(total_time: float, 
                   query_records: List[Tuple[float, int, List[Tuple[int, ...]]]]) -> None:
    """分析并显示测试结果"""
    # 计算统计信息
    query_times = [record[0] for record in query_records]
    avg_time = sum(query_times) / len(query_times)
    max_time = max(query_times)
    min_time = min(query_times)
    
    # 计算有结果的查询数量
    queries_with_results = sum(1 for _, _, results in query_records if results)
    
    print("\n====== 测试结果 ======")
    print(f"总查询数量: {len(query_records)}")
    print(f"总用时: {total_time:.2f} 秒")
    print(f"平均单次查询用时: {avg_time * 1000:.3f} 毫秒")
    print(f"最长单次查询用时: {max_time * 1000:.3f} 毫秒")
    print(f"最短单次查询用时: {min_time * 1000:.3f} 毫秒")
    print(f"每秒查询次数: {len(query_records) / total_time:.2f}")
    print(f"有结果的查询数量: {queries_with_results}")
    print(f"查询命中率: {(queries_with_results / len(query_records)) * 100:.2f}%")

def main():
    """主函数"""
    try:
        print("开始执行随机查询测试...")
        total_time, query_records = perform_random_queries()
        analyze_results(total_time, query_records)
        
    except Exception as e:
        print(f"程序执行出错: {str(e)}")

if __name__ == "__main__":
    main()
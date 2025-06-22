import sqlite3
import os
from typing import List, Tuple

def ensure_directory_exists(file_path: str) -> None:
    """确保目录存在，如果不存在则创建"""
    directory = os.path.dirname(file_path)
    if directory and not os.path.exists(directory):
        os.makedirs(directory)

def read_data_from_file(file_path: str) -> List[Tuple[int, ...]]:
    """从文件中读取数据"""
    data = []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                values = line.strip().split(':')
                if len(values) == 4:
                    data.append(tuple(map(int, values)))
    except FileNotFoundError:
        print(f"错误：找不到文件 {file_path}")
        raise
    except ValueError as e:
        print(f"错误：文件格式不正确，请确保所有值都是数字\n{str(e)}")
        raise
    return data

def create_database(input_file: str = "test/mock_interim.txt", 
                   db_file: str = "test/sqlite3.db") -> None:
    """创建SQLite数据库并建立表和索引"""
    ensure_directory_exists(db_file)
    
    if os.path.exists(db_file):
        os.remove(db_file)
    
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    
    # 创建表
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS data_table (
            col1 INTEGER,
            col2 INTEGER,
            col3 INTEGER,
            col4 INTEGER
        )
    ''')
    
    # 创建索引
    cursor.execute('CREATE INDEX idx_col2 ON data_table(col2)')
    
    try:
        # 读取并插入数据
        data = read_data_from_file(input_file)
        cursor.executemany('INSERT INTO data_table VALUES (?, ?, ?, ?)', data)
        conn.commit()
        
        # 验证数据
        cursor.execute('SELECT COUNT(*) FROM data_table')
        count = cursor.fetchone()[0]
        print(f"成功插入 {count} 条数据")
        
        print(f"数据库文件已创建: {os.path.abspath(db_file)}")
            
    except Exception as e:
        print(f"处理数据时出错: {str(e)}")
        conn.rollback()
        raise
    finally:
        cursor.close()
        conn.close()

if __name__ == "__main__":
    try:
        create_database()
        print("数据库创建成功，并已建立索引！")
    except Exception as e:
        print(f"程序执行出错: {str(e)}")
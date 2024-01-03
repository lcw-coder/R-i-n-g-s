def process_lammpstrj_file(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    frames = []
    current_frame = []

    # 分割每个帧并存储
    for line in lines:
        if line.startswith("ITEM: TIMESTEP"):
            if current_frame:
                frames.append(current_frame)
                current_frame = []
        current_frame.append(line)
    frames.append(current_frame)  # 添加最后一个帧

    # 处理每个帧的数据并写入一个 XYZ 格式文件
    with open(output_file, 'w') as xyz_file:
        for frame in frames:
            atom_data = frame[frame.index("ITEM: ATOMS element id q x y z type \n") + 1:]
            xyz_data = [f"{len(atom_data)}\nConverted from LAMMPS data"]
            for line in atom_data:
                columns = line.split()
                if len(columns) >= 7:
                    atom_xyz = ' '.join(columns[i] for i in [0, 3, 4, 5])
                    xyz_data.append(atom_xyz)
            xyz_file.write('\n'.join(xyz_data) + '\n')

# 处理第一个文件
process_lammpstrj_file('mydump1.lammpstrj', 'output1.xyz')

# 处理第二个文件
process_lammpstrj_file('mydump2.lammpstrj', 'output2.xyz')

# 读取第一个文件
with open('output1.xyz', 'r') as file1:
    content_file1 = file1.read()

# 读取第二个文件
with open('output2.xyz', 'r') as file2:
    content_file2 = file2.read()

# 合并内容
merged_content = content_file1 + content_file2

# 将合并后的内容写入新文件
with open('process.xyz', 'w') as merged_file:
    merged_file.write(merged_content)

##replace the O atoms into C in the xyz file
##writed by lcw
##2023-02-18
## Y11190003@mail.ecust.edu.cn

from tqdm import tqdm

# with open('output2.xyz', 'r') as f_in:
with open('process.xyz', 'r') as f_in:
    lines = f_in.readlines()

new_lines = []
for line in tqdm(lines, desc='Processing', unit=' line'):
    new_line = line.replace('O', 'C', 1).replace('N', 'C', 1)  # 替换第一个出现的'O'或'N'为'C'
    new_lines.append(new_line)

    # total_lines = sum(1 for _ in f_in) #获取总行数用于设置进度条
    # f_in.seek(0)  # 重置文件指针到文件开头
    # for line in tqdm(f_in, total=total_lines,desc="Processing", unit="line"):

    # new_line = line.replace('O', 'C', 1) # 只替换第一个C
    # new_line = line.replace('O', 'C', 1).replace('N', 'C', 1)  # 替换第一个出现的'O'或'N'为'C'
    # f_out.write(new_line)

with open('process2C.xyz', 'w') as f_out:
    f_out.writelines(new_lines)

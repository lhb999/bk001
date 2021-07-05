filename = "higgs-activity_time_1.txt"
old = open(filename, 'r')
tmp = []
data_list = []
line_count = 0
m_count = 0
edge = ("v ")
label = "0"

v_to_i = {}

def vertex_to_idx(v):
    if v not in v_to_i:
        v_to_i[v] = len(v_to_i)

lines = old.readlines()
for line in lines:
    line_count += 1
    data_list.append(line)
    spl1 = line.split(" ")

    vertex_to_idx(spl1[0])
    vertex_to_idx(spl1[1])

    tmp.append(spl1[0])
    tmp.append(spl1[1])

    if (line_count % 2000 == 0):
        distincted = sorted(list(set(tmp)))
        vertex_fileout = f"splited_{m_count}.txt"
        new = open(vertex_fileout, 'w')

        # 정점 먼저 쓰고
        for v in distincted:
            vv = v_to_i[v]
            new.write(f"v {vv} 0\n")

        # 엣지 쓰고
        for edge in data_list:
            spl = edge.split(" ")
            src = v_to_i[spl[0]]
            tgt = v_to_i[spl[1]]
            lbl = spl[2]
            str_with_e = f"e {src} {tgt} {lbl}"
            new.write(str_with_e)

        m_count += 1
        new.close()

        if line_count == 10000:
            break

old.close()

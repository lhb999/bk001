import numpy as np
import random

src = 1 # SRC vertex 초기값
tgt = 1000 # TGT vertex 초기값
size = 100 # 패턴 변화를 몇줄씩 줄지
line = 100
# 데이터셋을 몇줄 생성할건지
# 뭔진 모르겠지만 여기서 살짝 오류가있는게 1000* 1000의 라인이 출력 1000을 넣으면 10만줄,
# 100을 넣으면 만줄으로 출력.. 예측으로는 for문을 돌면서 size와 line의 갯수가 곱해지는듯 합니다!

"""돌아가는 순서
* SRC, TGT, LABLE의 크기는 같아야함
ex) size가 10이고 src는 1~10 tgt는 1000~1010의 범위를 가질때,
SRC, TGT는 범위의 값으로 랜덤 생성 (중복허용)
src = [1, 3, 4, 2, 5, 6, 10, 9, 8, 4]
tgt = [1001, 1008, 1008, 1005, 1006, 1000, 1009, 1009, 1007, 1002]
lable = [RT, MT, RE, RE, RT, RT, RT, RT, MT, RE] 이면

src[i] tgt[i] lable[i]
1, 1001, RT
3, 1008, MT
... 순으로 생성

"""
# lable 종류 생성
lable = ("RT", "MT", "RE")
lable_list = []

#lable을 랜덤으로 선택
for j in range(size):
    res = random.choice(lable)
    lable_list.append(res)
    # res = ''.join(lable_list)

#src와 tgt를 랜덤으로 생성
for i in range(line):
    with open('./rand.txt', 'a') as f:
        src_t = np.random.randint(src, src + 10, size)
        tgt_t = np.random.randint(tgt, tgt + 10, size)

        for x in range(0, len(tgt_t)):
            f.write(str(src_t[x]) + " " + str(tgt_t[x]) + " " + str(lable_list[x]) + "\n")

            if x == (size-1):
                src += 1
                tgt += 1
                continue

# line = 10000
#
# for i in range(line):
#     with open('./rand.txt','a') as f:
#         line_count = 0
#         src = 1
#         tgt = 1000
#
#         src_t = np.random.randint(src, src + 10, size=10)
#         tgt_t = np.random.randint(tgt, tgt + 10, size=10)
#         for x in range(0, len(tgt_t)):
#             f.write(str(src_t[x]) + " " + str(tgt_t[x]) + " " + str(lable_list[x]) + "\n")
# #
# # T = int(input())
# # 여러개의 테스트 케이스가 주어지므로, 각각을 처리합니다.
# line = "3 17 1 39 8 41 2 32 99 2"
#
#     n_sum = 0
#     nlist = line.split()
#     for v in nlist:
#         vv = int(v)
#         # print(vv)
#         if vv % 2 != 0:
#             n_sum += vv
#
#     print(f"#{1} {n_sum}")

"""
2
32850
160
"""
#
# hap = [32850, 160]
# don = [50000, 10000, 5000, 1000, 500, 100, 50, 10]
#
# for t in range(len(hap)):
#     tmp_hap = hap[t]
#     print(f"#{t+1}")
#     for d in don:
#         if tmp_hap >= d:
#             v = int(tmp_hap / d)
#             c = tmp_hap % d
#             tmp_hap = c
#             print(f"{v} ", end="")
#         else:
#             print(f"0 ", end="")
#     print()

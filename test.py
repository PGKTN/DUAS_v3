import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# CSV 파일에서 데이터 읽기
data = pd.read_csv("power_data.csv")

# 시간과 동력 데이터
time = data['Time']
power = data['Power']

# 각 구간 정의
intervals = []

# 시작 시간 및 동력 값 설정
start_time = time.iloc[0]
prev_power = power.iloc[0]
pattern_map = {}  # 동력 값에 따른 패턴 매핑을 위한 딕셔너리

# 각 데이터 포인트를 확인하면서 구간을 정의
for t, p in zip(time, power):
    if p != prev_power:  # 이전 동력 값과 현재 동력 값이 다르면 새로운 구간 시작
        intervals.append({
            "start": start_time,
            "end": t,
            "power": prev_power
        })
        start_time = t  # 새로운 구간의 시작 시간 설정
        prev_power = p  # 이전 동력 값 업데이트
    
    # 동력 값에 따른 패턴 매핑
    if p not in pattern_map:
        pattern_map[p] = len(pattern_map) + 1  # 동력 값에 따라 고유한 패턴 지정

# 마지막 구간 추가
intervals.append({
    "start": start_time,
    "end": time.iloc[-1],
    "power": prev_power
})

# 그래프 그리기
plt.figure(figsize=(10, 6))
ax = plt.gca()

# 각 구간에 대해 사각형 플롯
for interval in intervals:
    power = interval["power"]
    pattern = pattern_map[power]
    color = plt.cm.tab10(pattern)  # 동력 값에 따른 원색 사용
    rect = Rectangle((interval["start"], 0), interval["end"] - interval["start"], interval["power"], color=color, alpha=0.5, hatch='/' * pattern)
    ax.add_patch(rect)

# 그래프 제목과 라벨 설정
plt.title("Power Requirement over Time")
plt.xlabel("Time (seconds)")
plt.ylabel("Power (W)")

# x 축과 y 축 범위 설정
plt.xlim([time.min(), time.max()])
plt.ylim([0, max(power) * 1.1])

# 그래프 보이기
plt.grid(True)
plt.show()

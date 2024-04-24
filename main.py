import math

import numpy as np
import dem_data
from linked_line import LinkedLinePDE
from partition_algorithm import PartitionAlgorithm
from partition_optimal_reference_lineal_gorithm import PartitionOptimalRLAlgorithm

# main函数中声明的变量默认为global variable
# 而其他def函数中声明的变量则默认为local variable

# dem文件路径
file_path = r"C:\Users\GuoQing\Desktop\Experiment\S_PDERL\data\Copernicus_DSM_COG_10_N28_00_E097_00_DEM.tif"
# 分辨率
real_distance = 30
# 实例化dem类
dem = dem_data.Dem(file_path, real_distance)

# 保存对应象限需要添加的角度
quadrant = [0, math.pi * 3 / 2, math.pi, math.pi / 2]

# 结果矩阵 初始化为0
# 其中 result[x,y]表示的是第x行y列 对应在XY坐标系中是x->纵轴 y->横纵

# x_center = 97.4172
# y_center = 28.7107
# h_center = 3959.62948159039
# to_x = 97.418210
# to_y = 28.710065
# h_stand = 5

# radius = 128.84744002180742
# horizontal_start_angle = 0
# horizontal_end_angle = 360

"""
    pderl算法 
    并没有近似，对应的p和d均为观察点到各网格点的真实距离，其中x_grid_center仅用来对distance与height进行索引
"""


def analysis_by_pderl(x_center, y_center, h_center, to_x, to_y, h_stand):
    # 初始化
    dem_minx, dem_miny, per_x, per_y, x_grid_count, y_grid_count, see_height, x_grid_center, y_grid_center, min_x, max_x, min_y, max_y = get_initial_param(
        x_center, y_center, h_center, to_x, to_y, h_stand)
    x_distance, y_distance = get_xy_distance(x_center, y_center, dem_minx, dem_miny, x_grid_count, y_grid_count)

    # 纵向格网间距实地距离的倒数 用它来计算线段系数a
    u = 1 / (dem.dy * dem.rdy)

    # p,d根据不同区域会有所变化，但均符合以下两点：①p始终保持正值 ②d值随着调查线的构建应该逐渐变大
    # 邻近值p->1/x
    p = 0.0

    # 当前经纬度对应的高程值
    current_height = 0.0
    last_height = 0.0

    # 方向值d->y/x
    d = 0.0
    # 用于对比旧基准线上一次的当前点
    last_d = 0.0
    # 交点方向值
    cross_d = 0.0
    # 记录一次求交运算中交点区间的最小值
    min_d = 0.0

    # 线段对应系数a->u*(current_height-last_height)
    a = 0.0

    # 第一个点对应的e->height*p
    start_base_e = 0.0

    current_base_e = 0.0

    # 参考线的起始线段
    start_base_line = None
    # 用于对比旧基准线的当前点
    base_line = None
    # 即将接入新参考线的一部分线段,在不确定该线段定义域范围时暂存于该变量
    current_new_line = None

    # 格网点序号值,用来获取对应坐标的高程值height=dem[x_grid_index,y_grid_index]
    x_grid_index = 0
    y_grid_index = 0

    # x、y轴距离索引,用来获取对应地心坐标系x、y值 x=x_distance[x_distance_index]->将其转换为p、d、e
    x_distance_index = 0
    y_distance_index = 0

    # 边界索引值
    x_right_index = 2 * x_grid_count - 1
    y_top_index = 2 * y_grid_count - 1

    # 当前点对应参考线和调查线的e差值 dif_e=current_e(调查线)-current_base_e(参考线)
    e_diff = 0.0
    e_diff_last = 0.0

    # 右半面
    # 以下算法认为横向间隔和纵向间隔相同
    # 右半圆: 从南往北,从西往东算
    x_grid_index = x_grid_center + 1
    y_grid_index = min_y
    x_distance_index = x_grid_count
    y_distance_index = 0
    last_height = dem.height[x_grid_count, y_grid_count - 1] - see_height
    p = 1 / x_distance[x_distance_index]

    # 构造右半边初始参考线
    while y_grid_index <= max_y:
        # 当目标点的高度小于观察点的高度会存在负数的情况
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = y_distance[y_distance_index] * p
        a = u * (current_height - last_height)

        if start_base_line is None:
            # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
            start_base_line = LinkedLinePDE(d, a)
            base_line = start_base_line
            start_base_e = current_height * p
        else:
            base_line.link_forward(LinkedLinePDE(d, a))
            base_line = base_line.next

        #  第一列是可见的
        result[x_distance_index][y_distance_index] = 1
        last_height = current_height

        y_grid_index += 1
        y_distance_index += 1
    x_grid_index += 1
    x_distance_index += 1

    # 右边逐层计算参考线与点的显隐性
    while x_grid_index <= max_x:
        y_grid_index = min_y
        y_distance_index = 0

        # 分析区域外拓展了一点
        last_height = dem.height[x_grid_index, y_grid_index - 1] - see_height

        p = 1 / x_distance[x_distance_index]
        base_line = start_base_line

        # 判断起点的显隐性
        # 初始化起点的d、a
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = y_distance[y_distance_index] * p
        a = u * (current_height - last_height)

        """
        为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
        为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
        如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
        k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
        """
        while base_line.next.end_d < d:
            base_line = base_line.next
            # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
            start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        start_base_line = base_line
        current_base_e = start_base_e

        # 上一步相当于把延伸的点所参与的参考线段筛出去
        while base_line.end_d < d:
            base_line = base_line.next
            current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        """
        这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
        current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
        """
        e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

        # 可视性判断
        if e_diff >= 0:
            # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
            start_base_line = current_new_line = LinkedLinePDE(d, a)
            start_base_e = current_height * p
            result[x_distance_index][y_distance_index] = 1
        else:
            result[x_distance_index][y_distance_index] = 0
        last_d = d
        last_height = current_height

        # 从南向北过程
        y_grid_index += 1
        y_distance_index += 1

        # 判断后续点的显隐性别
        while y_grid_index <= max_y:
            # 形成新的调查线段
            d = y_distance[y_distance_index] * p
            current_height = dem.height[x_grid_index, y_grid_index] - see_height
            a = u * (current_height - last_height)

            min_d = last_d
            # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
            while base_line.end_d < d:
                """
                求出轴斜率差值的增量
                这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                """
                e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                # 判断是否为交点 如果前后可视性不同则比存在交点
                # 当前是可视的
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = base_line.end_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last
                min_d = base_line.end_d
                base_line = base_line.next

            e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
            if e_diff >= 0:
                # 变为不可视
                if e_diff_last < 0:
                    # d很小时算不准
                    if e_diff < 5e-15:
                        cross_d = d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                    # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                    current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                    current_new_line = current_new_line.next

                    # 这条语句同时还更新了参考线
                    current_new_line.link_forward(base_line)
            else:
                if e_diff_last >= 0:
                    if e_diff > -5e-15:
                        cross_d = min_d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    current_new_line = LinkedLinePDE(cross_d, base_line.a)
                    if base_line.pre is not None:
                        base_line.pre.link_forward(current_new_line)
                    else:
                        start_base_line = current_new_line
                        start_base_e = current_height * p - a * (d - cross_d)
            e_diff = e_diff_last

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                current_new_line.link_forward(LinkedLinePDE(d, a))
                current_new_line = current_new_line.next
                result[x_distance_index][y_distance_index] = 1
            else:
                result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从南向北过程
            y_grid_index += 1
            y_distance_index += 1

        # 最外层循环 对应网格点横坐标加1：从西向东过程
        x_grid_index += 1
        x_distance_index += 1
    # print("右面：\n", result)

    # 左半面
    # 左半圆:从北往南，从东往西算
    x_grid_index = x_grid_center
    y_grid_index = max_y
    x_distance_index = x_grid_count - 1
    y_distance_index = y_top_index
    last_height = dem.height[x_grid_index, y_grid_index + 1] - see_height
    start_base_line = None

    p = -1 / x_distance[x_distance_index]
    # 构造左半边初始参考线
    while y_grid_index >= min_y:
        # 当目标点的高度小于观察点的高度会存在负数的情况
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = -y_distance[y_distance_index] * p
        a = u * (current_height - last_height)

        if start_base_line is None:
            # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
            start_base_line = LinkedLinePDE(d, a)
            base_line = start_base_line
            start_base_e = current_height * p
        else:
            base_line.link_forward(LinkedLinePDE(d, a))
            base_line = base_line.next

        #  第一列是可见的
        result[x_distance_index][y_distance_index] = 1
        last_height = current_height

        y_grid_index -= 1
        y_distance_index -= 1

    x_grid_index -= 1
    x_distance_index -= 1

    # 左边逐层计算参考线与点的显隐性
    while x_grid_index >= min_x:
        y_grid_index = max_y
        y_distance_index = y_top_index

        # 分析区域外拓展了一点
        last_height = dem.height[x_grid_index, y_grid_index + 1] - see_height

        p = -1 / x_distance[x_distance_index]
        base_line = start_base_line

        # 判断起点的显隐性
        # 初始化起点的d、a
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = -y_distance[y_distance_index] * p
        a = u * (current_height - last_height)

        while base_line.next.end_d < d:
            base_line = base_line.next
            # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
            start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        start_base_line = base_line
        current_base_e = start_base_e

        # 上一步相当于把延伸的点所参与的参考线段筛出去
        while base_line.end_d < d:
            base_line = base_line.next
            current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        """
        这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
        current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
        """
        e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

        # 可视性判断
        if e_diff >= 0:
            # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
            start_base_line = current_new_line = LinkedLinePDE(d, a)
            start_base_e = current_height * p
            result[x_distance_index][y_distance_index] = 1
        else:
            result[x_distance_index][y_distance_index] = 0
        last_d = d
        last_height = current_height

        # 从北向南过程
        y_grid_index -= 1
        y_distance_index -= 1

        # 判断后续点的显隐性别
        while y_grid_index >= min_y:
            # 形成新的调查线段
            d = -y_distance[y_distance_index] * p
            current_height = dem.height[x_grid_index, y_grid_index] - see_height
            a = u * (current_height - last_height)

            min_d = last_d
            # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
            while base_line.end_d < d:
                """
                求出轴斜率差值的增量
                这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                """
                e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                # 判断是否为交点 如果前后可视性不同则比存在交点
                # 当前是可视的
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = base_line.end_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last
                min_d = base_line.end_d
                base_line = base_line.next

            e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
            if e_diff >= 0:
                # 变为不可视
                if e_diff_last < 0:
                    # d很小时算不准
                    if e_diff < 5e-15:
                        cross_d = d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                    # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                    current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                    current_new_line = current_new_line.next

                    # 这条语句同时还更新了参考线
                    current_new_line.link_forward(base_line)
            else:
                if e_diff_last >= 0:
                    if e_diff > -5e-15:
                        cross_d = min_d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    current_new_line = LinkedLinePDE(cross_d, base_line.a)
                    if base_line.pre is not None:
                        base_line.pre.link_forward(current_new_line)
                    else:
                        start_base_line = current_new_line
                        start_base_e = current_height * p - a * (d - cross_d)
            e_diff = e_diff_last

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                current_new_line.link_forward(LinkedLinePDE(d, a))
                current_new_line = current_new_line.next
                result[x_distance_index][y_distance_index] = 1
            else:
                result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从北向南过程
            y_grid_index -= 1
            y_distance_index -= 1

        # 最外层循环 对应网格点横坐标加1：从东向西过程
        x_grid_index -= 1
        x_distance_index -= 1
    # print("左面：\n", result)

    # 上半面及下半面->横向格网间距实地距离的倒数 用它来计算线段系数a
    u = 1 / (dem.dx * dem.rdx)

    # 上半面
    # 上半圆:从南往北，从东往西算
    x_grid_index = max_x
    y_grid_index = y_grid_center + 1
    x_distance_index = x_right_index
    y_distance_index = y_grid_count
    last_height = dem.height[x_grid_index + 1, y_grid_index] - see_height
    start_base_line = None

    p = 1 / y_distance[y_distance_index]

    # 构造上半边初始参考线
    while x_grid_index >= min_x:
        # 当目标点的高度小于观察点的高度会存在负数的情况
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = -x_distance[x_distance_index] * p
        a = u * (current_height - last_height)

        if start_base_line is None:
            # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
            start_base_line = LinkedLinePDE(d, a)
            base_line = start_base_line
            start_base_e = current_height * p
        else:
            base_line.link_forward(LinkedLinePDE(d, a))
            base_line = base_line.next

        #  第一列是可见的
        # result[x_distance_index][y_distance_index] = 1
        last_height = current_height

        x_grid_index -= 1
        x_distance_index -= 1
    y_grid_index += 1
    y_distance_index += 1

    # 上边逐层计算参考线与点的显隐性
    while y_grid_index <= max_y:
        x_grid_index = max_x
        x_distance_index = x_right_index

        # 分析区域外拓展了一点
        last_height = dem.height[x_grid_index + 1, y_grid_index] - see_height

        p = 1 / y_distance[y_distance_index]
        base_line = start_base_line

        # 判断起点的显隐性
        # 初始化起点的d、a
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = -x_distance[x_distance_index] * p
        a = u * (current_height - last_height)

        """
        为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
        为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
        如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
        k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
        """
        while base_line.next.end_d < d:
            base_line = base_line.next
            # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
            start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        start_base_line = base_line
        current_base_e = start_base_e

        # 上一步相当于把延伸的点所参与的参考线段筛出去
        while base_line.end_d < d:
            base_line = base_line.next
            current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        """
        这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
        current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
        """
        e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

        # 可视性判断
        if e_diff >= 0:
            # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
            start_base_line = current_new_line = LinkedLinePDE(d, a)
            start_base_e = current_height * p
            # result[x_distance_index][y_distance_index] = 1
        else:
            result[x_distance_index][y_distance_index] = 0
        last_d = d
        last_height = current_height

        # 从东向西过程
        x_grid_index -= 1
        x_distance_index -= 1

        # 判断后续点的显隐性别
        while x_grid_index >= min_x:
            # 形成新的调查线段
            d = -x_distance[x_distance_index] * p
            current_height = dem.height[x_grid_index, y_grid_index] - see_height
            a = u * (current_height - last_height)

            min_d = last_d
            # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
            while base_line.end_d < d:
                """
                求出轴斜率差值的增量
                这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                """
                e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                # 判断是否为交点 如果前后可视性不同则比存在交点
                # 当前是可视的
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = base_line.end_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last
                min_d = base_line.end_d
                base_line = base_line.next

            e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
            if e_diff >= 0:
                # 变为不可视
                if e_diff_last < 0:
                    # d很小时算不准
                    if e_diff < 5e-15:
                        cross_d = d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                    # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                    current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                    current_new_line = current_new_line.next

                    # 这条语句同时还更新了参考线
                    current_new_line.link_forward(base_line)
            else:
                if e_diff_last >= 0:
                    if e_diff > -5e-15:
                        cross_d = min_d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    current_new_line = LinkedLinePDE(cross_d, base_line.a)
                    if base_line.pre is not None:
                        base_line.pre.link_forward(current_new_line)
                    else:
                        start_base_line = current_new_line
                        start_base_e = current_height * p - a * (d - cross_d)
            e_diff = e_diff_last

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                current_new_line.link_forward(LinkedLinePDE(d, a))
                current_new_line = current_new_line.next

                # 上下半面 并没有赋值为1
                # result[x_distance_index][y_distance_index] = 1
            else:
                result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从东向西过程
            x_grid_index -= 1
            x_distance_index -= 1

        # 最外层循环 对应网格点横坐标加1：从南向北过程
        y_grid_index += 1
        y_distance_index += 1
    # print("上面：\n", result)

    # 下半面
    # 下半圆：从北往南，从西往东算
    x_grid_index = min_x
    y_grid_index = y_grid_center
    x_distance_index = 0
    y_distance_index = y_grid_count - 1
    last_height = dem.height[x_grid_index - 1, y_grid_index] - see_height
    start_base_line = None

    p = -1 / y_distance[y_distance_index]
    # 构造下半边初始参考线
    while x_grid_index <= max_x:
        # 当目标点的高度小于观察点的高度会存在负数的情况
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = x_distance[x_distance_index] * p
        a = u * (current_height - last_height)

        if start_base_line is None:
            # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
            start_base_line = LinkedLinePDE(d, a)
            base_line = start_base_line
            start_base_e = current_height * p
        else:
            base_line.link_forward(LinkedLinePDE(d, a))
            base_line = base_line.next

        #  第一列是可见的
        # result[x_distance_index][y_distance_index] = 1
        last_height = current_height

        x_grid_index += 1
        x_distance_index += 1

    y_grid_index -= 1
    y_distance_index -= 1

    # 下边逐层计算参考线与点的显隐性
    while y_grid_index >= min_y:
        x_grid_index = min_x
        x_distance_index = 0

        # 分析区域外拓展了一点
        last_height = dem.height[x_grid_index - 1, y_grid_index] - see_height

        p = -1 / y_distance[y_distance_index]
        base_line = start_base_line

        # 判断起点的显隐性
        # 初始化起点的d、a
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = x_distance[x_distance_index] * p
        a = u * (current_height - last_height)

        """
        为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
        为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
        如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
        k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
        """
        while base_line.next.end_d < d:
            base_line = base_line.next
            # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
            start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        start_base_line = base_line
        current_base_e = start_base_e

        # 上一步相当于把延伸的点所参与的参考线段筛出去
        while base_line.end_d < d:
            base_line = base_line.next
            current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        """
        这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
        current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
        """
        e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

        # 可视性判断
        if e_diff >= 0:
            # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
            start_base_line = current_new_line = LinkedLinePDE(d, a)
            start_base_e = current_height * p
            # result[x_distance_index][y_distance_index] = 1
        else:
            result[x_distance_index][y_distance_index] = 0

        last_d = d
        last_height = current_height

        # 从西向东过程
        x_grid_index += 1
        x_distance_index += 1

        # 判断后续点的显隐性别
        while x_grid_index <= max_x:
            # 形成新的调查线段
            d = x_distance[x_distance_index] * p
            current_height = dem.height[x_grid_index, y_grid_index] - see_height
            a = u * (current_height - last_height)

            min_d = last_d
            # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
            while base_line.end_d < d:
                """
                求出轴斜率差值的增量
                这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                """
                e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                # 判断是否为交点 如果前后可视性不同则比存在交点
                # 当前是可视的
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = base_line.end_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last
                min_d = base_line.end_d
                base_line = base_line.next

            e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
            if e_diff >= 0:
                # 变为不可视
                if e_diff_last < 0:
                    # d很小时算不准
                    if e_diff < 5e-15:
                        cross_d = d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                    # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                    current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                    current_new_line = current_new_line.next

                    # 这条语句同时还更新了参考线
                    current_new_line.link_forward(base_line)
            else:
                if e_diff_last >= 0:
                    if e_diff > -5e-15:
                        cross_d = min_d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    current_new_line = LinkedLinePDE(cross_d, base_line.a)
                    if base_line.pre is not None:
                        base_line.pre.link_forward(current_new_line)
                    else:
                        start_base_line = current_new_line
                        start_base_e = current_height * p - a * (d - cross_d)
            e_diff = e_diff_last

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                current_new_line.link_forward(LinkedLinePDE(d, a))
                current_new_line = current_new_line.next

                # result[x_distance_index][y_distance_index] = 1
            else:
                result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从西向东过程
            x_grid_index += 1
            x_distance_index += 1

        # 最外层循环 对应网格点横坐标加1：从北向南过程
        y_grid_index -= 1
        y_distance_index -= 1

    # print("最终结果:\n", result)

def analysis_by_pderl_with_angle(x_center, y_center, h_center, radius, h_stand, horizontal_start_angle,
                                 horizontal_end_angle, current_dem=dem):
    # 初始化
    global dem
    dem = current_dem
    # 初始化
    x_grid_observe, y_grid_observe, dem_minx, dem_miny, x_grid_count, y_grid_count, see_height, x_grid_center, y_grid_center, min_x, max_x, min_y, max_y, is_in_grid_point, result = get_initial_param_spderl(
        x_center, y_center, h_center, radius, h_stand)
    x_distance, y_distance = get_xy_distance(x_center, y_center, dem_minx, dem_miny, x_grid_count, y_grid_count)

    # 纵向格网间距实地距离的倒数 用它来计算线段系数a
    u = 1 / (dem.dy * dem.rdy)

    # p,d根据不同区域会有所变化，但均符合以下两点：①p始终保持正值 ②d值随着调查线的构建应该逐渐变大
    # 邻近值p->1/x
    p = 0.0

    # 当前经纬度对应的高程值
    current_height = 0.0
    last_height = 0.0

    # 方向值d->y/x
    d = 0.0
    # 用于对比旧基准线上一次的当前点
    last_d = 0.0
    # 交点方向值
    cross_d = 0.0
    # 记录一次求交运算中交点区间的最小值
    min_d = 0.0

    # 线段对应系数a->u*(current_height-last_height)
    a = 0.0

    # 第一个点对应的e->height*p
    start_base_e = 0.0

    current_base_e = 0.0

    # 参考线的起始线段
    start_base_line = None
    # 用于对比旧基准线的当前点
    base_line = None
    # 即将接入新参考线的一部分线段,在不确定该线段定义域范围时暂存于该变量
    current_new_line = None

    # 格网点序号值,用来获取对应坐标的高程值height=dem[x_grid_index,y_grid_index]
    x_grid_index = 0
    y_grid_index = 0

    # x、y轴距离索引,用来获取对应地心坐标系x、y值 x=x_distance[x_distance_index]->将其转换为p、d、e
    x_distance_index = 0
    y_distance_index = 0

    # 边界索引值
    x_right_index = 2 * x_grid_count - 1
    y_top_index = 2 * y_grid_count - 1

    # 当前点对应参考线和调查线的e差值 dif_e=current_e(调查线)-current_base_e(参考线)
    e_diff = 0.0
    e_diff_last = 0.0

    # 右半面
    # 以下算法认为横向间隔和纵向间隔相同
    # 右半圆: 从南往北,从西往东算
    x_grid_index = x_grid_center + 1
    y_grid_index = min_y
    x_distance_index = x_grid_count
    y_distance_index = 0
    last_height = dem.height[x_grid_count, y_grid_count - 1] - see_height
    p = 1 / x_distance[x_distance_index]

    # 构造右半边初始参考线
    while y_grid_index <= max_y:
        # 当目标点的高度小于观察点的高度会存在负数的情况
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = y_distance[y_distance_index] * p
        a = u * (current_height - last_height)

        if start_base_line is None:
            # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
            start_base_line = LinkedLinePDE(d, a)
            base_line = start_base_line
            start_base_e = current_height * p
        else:
            base_line.link_forward(LinkedLinePDE(d, a))
            base_line = base_line.next

        #  第一列是可见的
        result[x_distance_index][y_distance_index] = 1
        last_height = current_height

        y_grid_index += 1
        y_distance_index += 1
    x_grid_index += 1
    x_distance_index += 1

    # 右边逐层计算参考线与点的显隐性
    while x_grid_index <= max_x:
        y_grid_index = min_y
        y_distance_index = 0

        # 分析区域外拓展了一点
        last_height = dem.height[x_grid_index, y_grid_index - 1] - see_height

        p = 1 / x_distance[x_distance_index]
        base_line = start_base_line

        # 判断起点的显隐性
        # 初始化起点的d、a
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = y_distance[y_distance_index] * p
        a = u * (current_height - last_height)

        """
        为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
        为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
        如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
        k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
        """
        while base_line.next.end_d < d:
            base_line = base_line.next
            # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
            start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        start_base_line = base_line
        current_base_e = start_base_e

        # 上一步相当于把延伸的点所参与的参考线段筛出去
        while base_line.end_d < d:
            base_line = base_line.next
            current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        """
        这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
        current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
        """
        e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

        # 可视性判断
        if e_diff >= 0:
            # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
            start_base_line = current_new_line = LinkedLinePDE(d, a)
            start_base_e = current_height * p
            result[x_distance_index][y_distance_index] = 1
        else:
            result[x_distance_index][y_distance_index] = 0
        last_d = d
        last_height = current_height

        # 从南向北过程
        y_grid_index += 1
        y_distance_index += 1

        # 判断后续点的显隐性别
        while y_grid_index <= max_y:
            # 形成新的调查线段
            d = y_distance[y_distance_index] * p
            current_height = dem.height[x_grid_index, y_grid_index] - see_height
            a = u * (current_height - last_height)

            min_d = last_d
            # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
            while base_line.end_d < d:
                """
                求出轴斜率差值的增量
                这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                """
                e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                # 判断是否为交点 如果前后可视性不同则比存在交点
                # 当前是可视的
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = base_line.end_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last
                min_d = base_line.end_d
                base_line = base_line.next

            e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
            if e_diff >= 0:
                # 变为不可视
                if e_diff_last < 0:
                    # d很小时算不准
                    if e_diff < 5e-15:
                        cross_d = d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                    # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                    current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                    current_new_line = current_new_line.next

                    # 这条语句同时还更新了参考线
                    current_new_line.link_forward(base_line)
            else:
                if e_diff_last >= 0:
                    if e_diff > -5e-15:
                        cross_d = min_d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    current_new_line = LinkedLinePDE(cross_d, base_line.a)
                    if base_line.pre is not None:
                        base_line.pre.link_forward(current_new_line)
                    else:
                        start_base_line = current_new_line
                        start_base_e = current_height * p - a * (d - cross_d)
            e_diff = e_diff_last

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                current_new_line.link_forward(LinkedLinePDE(d, a))
                current_new_line = current_new_line.next
                current_quadrant = 1 if y_grid_index >= y_grid_center else 4
                is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                                y_grid_observe, x_grid_index, y_grid_index, current_quadrant)
                if is_in_range and x_distance_index >= x_grid_count:
                    result[x_distance_index][y_distance_index] = 1
            else:
                result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从南向北过程
            y_grid_index += 1
            y_distance_index += 1

        # 最外层循环 对应网格点横坐标加1：从西向东过程
        x_grid_index += 1
        x_distance_index += 1
    # print("右面：\n", result)

    # 左半面
    # 左半圆:从北往南，从东往西算
    x_grid_index = x_grid_center
    y_grid_index = max_y
    x_distance_index = x_grid_count - 1
    y_distance_index = y_top_index
    last_height = dem.height[x_grid_index, y_grid_index + 1] - see_height
    start_base_line = None

    p = -1 / x_distance[x_distance_index]
    # 构造左半边初始参考线
    while y_grid_index >= min_y:
        # 当目标点的高度小于观察点的高度会存在负数的情况
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = -y_distance[y_distance_index] * p
        a = u * (current_height - last_height)

        if start_base_line is None:
            # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
            start_base_line = LinkedLinePDE(d, a)
            base_line = start_base_line
            start_base_e = current_height * p
        else:
            base_line.link_forward(LinkedLinePDE(d, a))
            base_line = base_line.next

        #  第一列是可见的
        result[x_distance_index][y_distance_index] = 1
        last_height = current_height

        y_grid_index -= 1
        y_distance_index -= 1

    x_grid_index -= 1
    x_distance_index -= 1

    # 左边逐层计算参考线与点的显隐性
    while x_grid_index >= min_x:
        y_grid_index = max_y
        y_distance_index = y_top_index

        # 分析区域外拓展了一点
        last_height = dem.height[x_grid_index, y_grid_index + 1] - see_height

        p = -1 / x_distance[x_distance_index]
        base_line = start_base_line

        # 判断起点的显隐性
        # 初始化起点的d、a
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = -y_distance[y_distance_index] * p
        a = u * (current_height - last_height)

        while base_line.next.end_d < d:
            base_line = base_line.next
            # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
            start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        start_base_line = base_line
        current_base_e = start_base_e

        # 上一步相当于把延伸的点所参与的参考线段筛出去
        while base_line.end_d < d:
            base_line = base_line.next
            current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        """
        这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
        current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
        """
        e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

        # 可视性判断
        if e_diff >= 0:
            # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
            start_base_line = current_new_line = LinkedLinePDE(d, a)
            start_base_e = current_height * p
            result[x_distance_index][y_distance_index] = 1
        else:
            result[x_distance_index][y_distance_index] = 0
        last_d = d
        last_height = current_height

        # 从北向南过程
        y_grid_index -= 1
        y_distance_index -= 1

        # 判断后续点的显隐性别
        while y_grid_index >= min_y:
            # 形成新的调查线段
            d = -y_distance[y_distance_index] * p
            current_height = dem.height[x_grid_index, y_grid_index] - see_height
            a = u * (current_height - last_height)

            min_d = last_d
            # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
            while base_line.end_d < d:
                """
                求出轴斜率差值的增量
                这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                """
                e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                # 判断是否为交点 如果前后可视性不同则比存在交点
                # 当前是可视的
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = base_line.end_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last
                min_d = base_line.end_d
                base_line = base_line.next

            e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
            if e_diff >= 0:
                # 变为不可视
                if e_diff_last < 0:
                    # d很小时算不准
                    if e_diff < 5e-15:
                        cross_d = d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                    # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                    current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                    current_new_line = current_new_line.next

                    # 这条语句同时还更新了参考线
                    current_new_line.link_forward(base_line)
            else:
                if e_diff_last >= 0:
                    if e_diff > -5e-15:
                        cross_d = min_d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    current_new_line = LinkedLinePDE(cross_d, base_line.a)
                    if base_line.pre is not None:
                        base_line.pre.link_forward(current_new_line)
                    else:
                        start_base_line = current_new_line
                        start_base_e = current_height * p - a * (d - cross_d)
            e_diff = e_diff_last

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                current_new_line.link_forward(LinkedLinePDE(d, a))
                current_new_line = current_new_line.next
                current_quadrant = 2 if y_grid_index >= y_grid_center else 3
                is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                                y_grid_observe, x_grid_index, y_grid_index, current_quadrant)
                if is_in_range and x_distance_index <= x_grid_count:
                    result[x_distance_index][y_distance_index] = 1
            else:
                result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从北向南过程
            y_grid_index -= 1
            y_distance_index -= 1

        # 最外层循环 对应网格点横坐标加1：从东向西过程
        x_grid_index -= 1
        x_distance_index -= 1
    # print("左面：\n", result)

    # 上半面及下半面->横向格网间距实地距离的倒数 用它来计算线段系数a
    u = 1 / (dem.dx * dem.rdx)

    # 上半面
    # 上半圆:从南往北，从东往西算
    x_grid_index = max_x
    y_grid_index = y_grid_center + 1
    x_distance_index = x_right_index
    y_distance_index = y_grid_count
    last_height = dem.height[x_grid_index + 1, y_grid_index] - see_height
    start_base_line = None

    p = 1 / y_distance[y_distance_index]

    # 构造上半边初始参考线
    while x_grid_index >= min_x:
        # 当目标点的高度小于观察点的高度会存在负数的情况
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = -x_distance[x_distance_index] * p
        a = u * (current_height - last_height)

        if start_base_line is None:
            # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
            start_base_line = LinkedLinePDE(d, a)
            base_line = start_base_line
            start_base_e = current_height * p
        else:
            base_line.link_forward(LinkedLinePDE(d, a))
            base_line = base_line.next

        #  第一列是可见的
        # result[x_distance_index][y_distance_index] = 1
        last_height = current_height

        x_grid_index -= 1
        x_distance_index -= 1
    y_grid_index += 1
    y_distance_index += 1

    # 上边逐层计算参考线与点的显隐性
    while y_grid_index <= max_y:
        x_grid_index = max_x
        x_distance_index = x_right_index

        # 分析区域外拓展了一点
        last_height = dem.height[x_grid_index + 1, y_grid_index] - see_height

        p = 1 / y_distance[y_distance_index]
        base_line = start_base_line

        # 判断起点的显隐性
        # 初始化起点的d、a
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = -x_distance[x_distance_index] * p
        a = u * (current_height - last_height)

        """
        为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
        为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
        如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
        k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
        """
        while base_line.next.end_d < d:
            base_line = base_line.next
            # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
            start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        start_base_line = base_line
        current_base_e = start_base_e

        # 上一步相当于把延伸的点所参与的参考线段筛出去
        while base_line.end_d < d:
            base_line = base_line.next
            current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        """
        这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
        current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
        """
        e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

        # 可视性判断
        if e_diff >= 0:
            # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
            start_base_line = current_new_line = LinkedLinePDE(d, a)
            start_base_e = current_height * p
            # result[x_distance_index][y_distance_index] = 1
        else:
            result[x_distance_index][y_distance_index] = 0
        last_d = d
        last_height = current_height

        # 从东向西过程
        x_grid_index -= 1
        x_distance_index -= 1

        # 判断后续点的显隐性别
        while x_grid_index >= min_x:
            # 形成新的调查线段
            d = -x_distance[x_distance_index] * p
            current_height = dem.height[x_grid_index, y_grid_index] - see_height
            a = u * (current_height - last_height)

            min_d = last_d
            # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
            while base_line.end_d < d:
                """
                求出轴斜率差值的增量
                这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                """
                e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                # 判断是否为交点 如果前后可视性不同则比存在交点
                # 当前是可视的
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = base_line.end_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last
                min_d = base_line.end_d
                base_line = base_line.next

            e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
            if e_diff >= 0:
                # 变为不可视
                if e_diff_last < 0:
                    # d很小时算不准
                    if e_diff < 5e-15:
                        cross_d = d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                    # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                    current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                    current_new_line = current_new_line.next

                    # 这条语句同时还更新了参考线
                    current_new_line.link_forward(base_line)
            else:
                if e_diff_last >= 0:
                    if e_diff > -5e-15:
                        cross_d = min_d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    current_new_line = LinkedLinePDE(cross_d, base_line.a)
                    if base_line.pre is not None:
                        base_line.pre.link_forward(current_new_line)
                    else:
                        start_base_line = current_new_line
                        start_base_e = current_height * p - a * (d - cross_d)
            e_diff = e_diff_last

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                current_new_line.link_forward(LinkedLinePDE(d, a))
                current_new_line = current_new_line.next

                # 上下半面 并没有赋值为1
                # result[x_distance_index][y_distance_index] = 1
            else:
                result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从东向西过程
            x_grid_index -= 1
            x_distance_index -= 1

        # 最外层循环 对应网格点横坐标加1：从南向北过程
        y_grid_index += 1
        y_distance_index += 1
    # print("上面：\n", result)

    # 下半面
    # 下半圆：从北往南，从西往东算
    x_grid_index = min_x
    y_grid_index = y_grid_center
    x_distance_index = 0
    y_distance_index = y_grid_count - 1
    last_height = dem.height[x_grid_index - 1, y_grid_index] - see_height
    start_base_line = None

    p = -1 / y_distance[y_distance_index]
    # 构造下半边初始参考线
    while x_grid_index <= max_x:
        # 当目标点的高度小于观察点的高度会存在负数的情况
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = x_distance[x_distance_index] * p
        a = u * (current_height - last_height)

        if start_base_line is None:
            # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
            start_base_line = LinkedLinePDE(d, a)
            base_line = start_base_line
            start_base_e = current_height * p
        else:
            base_line.link_forward(LinkedLinePDE(d, a))
            base_line = base_line.next

        #  第一列是可见的
        # result[x_distance_index][y_distance_index] = 1
        last_height = current_height

        x_grid_index += 1
        x_distance_index += 1

    y_grid_index -= 1
    y_distance_index -= 1

    # 下边逐层计算参考线与点的显隐性
    while y_grid_index >= min_y:
        x_grid_index = min_x
        x_distance_index = 0

        # 分析区域外拓展了一点
        last_height = dem.height[x_grid_index - 1, y_grid_index] - see_height

        p = -1 / y_distance[y_distance_index]
        base_line = start_base_line

        # 判断起点的显隐性
        # 初始化起点的d、a
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = x_distance[x_distance_index] * p
        a = u * (current_height - last_height)

        """
        为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
        为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
        如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
        k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
        """
        while base_line.next.end_d < d:
            base_line = base_line.next
            # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
            start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        start_base_line = base_line
        current_base_e = start_base_e

        # 上一步相当于把延伸的点所参与的参考线段筛出去
        while base_line.end_d < d:
            base_line = base_line.next
            current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        """
        这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
        current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
        """
        e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

        # 可视性判断
        if e_diff >= 0:
            # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
            start_base_line = current_new_line = LinkedLinePDE(d, a)
            start_base_e = current_height * p
            # result[x_distance_index][y_distance_index] = 1
        else:
            result[x_distance_index][y_distance_index] = 0

        last_d = d
        last_height = current_height

        # 从西向东过程
        x_grid_index += 1
        x_distance_index += 1

        # 判断后续点的显隐性别
        while x_grid_index <= max_x:
            # 形成新的调查线段
            d = x_distance[x_distance_index] * p
            current_height = dem.height[x_grid_index, y_grid_index] - see_height
            a = u * (current_height - last_height)

            min_d = last_d
            # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
            while base_line.end_d < d:
                """
                求出轴斜率差值的增量
                这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                """
                e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                # 判断是否为交点 如果前后可视性不同则比存在交点
                # 当前是可视的
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = base_line.end_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last
                min_d = base_line.end_d
                base_line = base_line.next

            e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
            if e_diff >= 0:
                # 变为不可视
                if e_diff_last < 0:
                    # d很小时算不准
                    if e_diff < 5e-15:
                        cross_d = d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                    # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                    current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                    current_new_line = current_new_line.next

                    # 这条语句同时还更新了参考线
                    current_new_line.link_forward(base_line)
            else:
                if e_diff_last >= 0:
                    if e_diff > -5e-15:
                        cross_d = min_d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    current_new_line = LinkedLinePDE(cross_d, base_line.a)
                    if base_line.pre is not None:
                        base_line.pre.link_forward(current_new_line)
                    else:
                        start_base_line = current_new_line
                        start_base_e = current_height * p - a * (d - cross_d)
            e_diff = e_diff_last

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                current_new_line.link_forward(LinkedLinePDE(d, a))
                current_new_line = current_new_line.next

                # result[x_distance_index][y_distance_index] = 1
            else:
                result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从西向东过程
            x_grid_index += 1
            x_distance_index += 1

        # 最外层循环 对应网格点横坐标加1：从北向南过程
        y_grid_index -= 1
        y_distance_index -= 1

    return result
    # print("最终结果:\n", result)

    # 比较与原算法结果是否有异
    true_result = judge_with_true.get_true_result()
    # print("是否相同：\n", judge_with_true.are_arrays_equal(result, true_result))


def analysis_by_xpderl(x_center, y_center, h_center, to_x, to_y, h_stand):
    # 初始化
    dem_minx, dem_miny, per_x, per_y, x_grid_count, y_grid_count, see_height, x_grid_center, y_grid_center, min_x, max_x, min_y, max_y = get_initial_param(
        x_center, y_center, h_center, to_x, to_y, h_stand)
    x_distance, y_distance = get_xy_distance(x_center, y_center, dem_minx, dem_miny, x_grid_count, y_grid_count)

    # 纵向格网间距实地距离的倒数 用它来计算线段系数a
    u = 1 / (dem.dy * dem.rdy)

    # p,d根据不同区域会有所变化，但均符合以下两点：①p始终保持正值 ②d值随着调查线的构建应该逐渐变大
    # 邻近值p->1/x
    p = 0.0

    # 当前经纬度对应的高程值
    current_height = 0.0
    last_height = 0.0

    # 方向值d->y/x
    d = 0.0
    # 用于对比旧基准线上一次的当前点
    last_d = 0.0
    # 交点方向值
    cross_d = 0.0
    # 记录一次求交运算中交点区间的最小值
    min_d = 0.0

    # 线段对应系数a->u*(current_height-last_height)
    a = 0.0

    # 第一个点对应的e->height*p
    start_base_e = 0.0

    current_base_e = 0.0

    # 参考线的起始线段
    start_base_line = None
    # 用于对比旧基准线的当前点
    base_line = None
    # 即将接入新参考线的一部分线段,在不确定该线段定义域范围时暂存于该变量
    current_new_line = None

    # 求出中心点在中心格网的相对位置 这里的中心网格不代表为整个结果的中心点 而是代表对一个网格框(由四个点组成)他们的中心
    # dem_min对应为一个网格点的经纬度

    lon_grid_center = (x_center - dem_minx) % dem.dx
    lat_grid_center = (y_center - dem_miny) % dem.dy

    # 针对不同区域来判断
    # 观察点处于对角线上==>分割线与 DEM 网格中的对角线重合，则无需扩展计算区域
    # 如果与某一条对角线重合也不用扩展
    # 拓展边界应进行调整 若在对角线上则无需拓展 反之拓展
    is_center_left_top = lat_grid_center / lon_grid_center - dem.dy / dem.dx >= 0
    # 结合图像的几何意义更易于理解
    is_center_left_bottom = dem.dy * lon_grid_center + dem.dx * lat_grid_center <= dem.dx * dem.dy

    # 格网点序号值,用来获取对应坐标的高程值height=dem[x_grid_index,y_grid_index]
    x_grid_index = 0
    y_grid_index = 0

    # x、y轴距离索引,用来获取对应地心坐标系x、y值 x=x_distance[x_distance_index]->将其转换为p、d、e
    x_distance_index = 0
    y_distance_index = 0

    # 边界索引值
    x_right_index = 2 * x_grid_count - 1
    y_top_index = 2 * y_grid_count - 1

    # 当前点对应参考线和调查线的e差值 dif_e=current_e(调查线)-current_base_e(参考线)
    e_diff = 0.0
    e_diff_last = 0.0

    # 右半面
    # 以下算法认为横向间隔和纵向间隔相同
    # 右半圆: 从南往北,从西往东算
    # 通过约束纵向y索引来限制判断区域
    # 拓展边界
    start_index_adjust = -1 if is_center_left_bottom else 0
    end_index_adjust = 1 if is_center_left_top else 0

    # 约束循环条件
    start_y_index = y_grid_center + start_index_adjust
    end_y_index = y_grid_center + end_index_adjust + 1

    # 调整对应的d边界
    start_distance_y_index = y_grid_count + start_index_adjust - 1

    x_grid_index = x_grid_center + 1
    y_grid_index = start_y_index
    x_distance_index = x_grid_count
    y_distance_index = start_distance_y_index
    last_height = dem.height[x_grid_index, y_grid_index - 1] - see_height
    p = 1 / x_distance[x_distance_index]

    # 观察点不是网格点该算法也只是将其视为中心网格点上方或者下方的网格点 并没有对其求值
    # 构造右半边初始参考线
    while y_grid_index <= end_y_index:
        # 当目标点的高度小于观察点的高度会存在负数的情况
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = y_distance[y_distance_index] * p
        a = u * (current_height - last_height)

        if start_base_line is None:
            # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
            start_base_line = LinkedLinePDE(d, a)
            base_line = start_base_line
            start_base_e = current_height * p
        else:
            base_line.link_forward(LinkedLinePDE(d, a))
            base_line = base_line.next

        #  第一列是可见的
        result[x_distance_index][y_distance_index] = 1
        last_height = current_height

        y_grid_index += 1
        y_distance_index += 1

    # 区域内第一个观察点的y坐标
    start_distance_y_index -= 1

    # 向y方向上下拓展边界 x型分区
    start_y_index -= 1
    end_y_index += 1

    # 目标点继续向右移动
    x_grid_index += 1
    x_distance_index += 1

    # 右边逐层计算参考线与点的显隐性
    while x_grid_index <= max_x:
        # 防止超出结果数组边界,回调索引
        if start_distance_y_index < 0:
            start_y_index += 1
            start_distance_y_index = 0

        y_grid_index = start_y_index
        y_distance_index = start_distance_y_index

        # 分析区域外拓展了一点
        last_height = dem.height[x_grid_index, y_grid_index - 1] - see_height

        p = 1 / x_distance[x_distance_index]
        base_line = start_base_line

        # 判断起点的显隐性
        # 初始化起点的d、a
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = y_distance[y_distance_index] * p
        a = u * (current_height - last_height)

        """
        为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
        为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
        如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
        k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
        """
        while base_line.next.end_d < d:
            base_line = base_line.next
            # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
            start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        start_base_line = base_line
        current_base_e = start_base_e

        # 上一步相当于把延伸的点所参与的参考线段筛出去
        while base_line.end_d < d:
            base_line = base_line.next
            current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        """
        这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
        current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
        """
        e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

        # 可视性判断
        if e_diff >= 0:
            # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
            start_base_line = current_new_line = LinkedLinePDE(d, a)
            start_base_e = current_height * p
            result[x_distance_index][y_distance_index] = 1
        else:
            result[x_distance_index][y_distance_index] = 0
        last_d = d
        last_height = current_height

        # 从南向北过程
        y_grid_index += 1
        y_distance_index += 1

        # 判断后续点的显隐性别
        while y_grid_index <= max_y and y_grid_index <= end_y_index:
            # 形成新的调查线段
            d = y_distance[y_distance_index] * p
            current_height = dem.height[x_grid_index, y_grid_index] - see_height
            a = u * (current_height - last_height)

            min_d = last_d
            # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
            while base_line.end_d < d:
                """
                求出轴斜率差值的增量
                这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                """
                e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                # 判断是否为交点 如果前后可视性不同则比存在交点
                # 当前是可视的
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = base_line.end_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last
                min_d = base_line.end_d
                base_line = base_line.next

            e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
            if e_diff >= 0:
                # 变为不可视
                if e_diff_last < 0:
                    # d很小时算不准
                    if e_diff < 5e-15:
                        cross_d = d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                    # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                    current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                    current_new_line = current_new_line.next

                    # 这条语句同时还更新了参考线
                    current_new_line.link_forward(base_line)
            else:
                if e_diff_last >= 0:
                    if e_diff > -5e-15:
                        cross_d = min_d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    current_new_line = LinkedLinePDE(cross_d, base_line.a)
                    if base_line.pre is not None:
                        base_line.pre.link_forward(current_new_line)
                    else:
                        start_base_line = current_new_line
                        start_base_e = current_height * p - a * (d - cross_d)
            e_diff = e_diff_last

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                current_new_line.link_forward(LinkedLinePDE(d, a))
                current_new_line = current_new_line.next
                result[x_distance_index][y_distance_index] = 1
            else:
                result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从南向北过程
            y_grid_index += 1
            y_distance_index += 1

        # 最外层循环 对应网格点横坐标加1：从西向东过程
        x_grid_index += 1
        x_distance_index += 1
        start_y_index -= 1
        start_distance_y_index -= 1
        end_y_index += 1
    # print("右面：\n", result)

    # 上半面及下半面->横向格网间距实地距离的倒数 用它来计算线段系数a
    u = 1 / (dem.dx * dem.rdx)

    # 上半面
    # 上半圆:从南往北，从东往西算
    is_center_right_bottom = lat_grid_center / lon_grid_center - dem.dy / dem.dx <= 0

    start_index_adjust = 1 if is_center_right_bottom else 0
    end_index_adjust = -1 if is_center_left_bottom else 0

    start_x_index = x_grid_center + start_index_adjust + 1
    start_distance_x_index = x_grid_count + start_index_adjust
    end_x_index = x_grid_center + end_index_adjust

    x_grid_index = start_x_index
    y_grid_index = y_grid_center + 1
    x_distance_index = start_distance_x_index
    y_distance_index = y_grid_count
    last_height = dem.height[x_grid_index + 1, y_grid_index] - see_height
    start_base_line = None

    p = 1 / y_distance[y_distance_index]

    # 构造上半边初始参考线
    while x_grid_index >= end_x_index:
        # 当目标点的高度小于观察点的高度会存在负数的情况
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = -x_distance[x_distance_index] * p
        a = u * (current_height - last_height)

        if start_base_line is None:
            # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
            start_base_line = LinkedLinePDE(d, a)
            base_line = start_base_line
            start_base_e = current_height * p
        else:
            base_line.link_forward(LinkedLinePDE(d, a))
            base_line = base_line.next

        #  第一列是可见的
        # result[x_distance_index][y_distance_index] = 1
        last_height = current_height

        x_grid_index -= 1
        x_distance_index -= 1

    # 拓展 x 方向
    start_x_index += 1
    start_distance_x_index += 1
    end_x_index -= 1

    y_grid_index += 1
    y_distance_index += 1

    # 上边逐层计算参考线与点的显隐性
    # 防止超出x边界
    max_right_index = 2 * x_grid_count - 1
    while y_grid_index <= max_y:
        if start_distance_x_index > max_right_index:
            start_x_index -= 1
            start_distance_x_index = max_right_index
        x_grid_index = start_x_index
        x_distance_index = start_distance_x_index

        # 分析区域外拓展了一点
        last_height = dem.height[x_grid_index + 1, y_grid_index] - see_height

        p = 1 / y_distance[y_distance_index]
        base_line = start_base_line

        # 判断起点的显隐性
        # 初始化起点的d、a
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = -x_distance[x_distance_index] * p
        a = u * (current_height - last_height)

        """
        为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
        为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
        如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
        k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
        """
        while base_line.next.end_d < d:
            base_line = base_line.next
            # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
            start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        start_base_line = base_line
        current_base_e = start_base_e

        # 上一步相当于把延伸的点所参与的参考线段筛出去
        while base_line.end_d < d:
            base_line = base_line.next
            current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        """
        这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
        current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
        """
        e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

        # 可视性判断
        if e_diff >= 0:
            # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
            start_base_line = current_new_line = LinkedLinePDE(d, a)
            start_base_e = current_height * p
            # result[x_distance_index][y_distance_index] = 1
        else:
            result[x_distance_index][y_distance_index] = 0
        last_d = d
        last_height = current_height

        # 从东向西过程
        x_grid_index -= 1
        x_distance_index -= 1

        # 判断后续点的显隐性别
        while x_grid_index >= min_x and x_grid_index >= end_x_index:
            # 形成新的调查线段
            d = -x_distance[x_distance_index] * p
            current_height = dem.height[x_grid_index, y_grid_index] - see_height
            a = u * (current_height - last_height)

            min_d = last_d
            # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
            while base_line.end_d < d:
                """
                求出轴斜率差值的增量
                这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                """
                e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                # 判断是否为交点 如果前后可视性不同则比存在交点
                # 当前是可视的
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = base_line.end_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last
                min_d = base_line.end_d
                base_line = base_line.next

            e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
            if e_diff >= 0:
                # 变为不可视
                if e_diff_last < 0:
                    # d很小时算不准
                    if e_diff < 5e-15:
                        cross_d = d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                    # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                    current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                    current_new_line = current_new_line.next

                    # 这条语句同时还更新了参考线
                    current_new_line.link_forward(base_line)
            else:
                if e_diff_last >= 0:
                    if e_diff > -5e-15:
                        cross_d = min_d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    current_new_line = LinkedLinePDE(cross_d, base_line.a)
                    if base_line.pre is not None:
                        base_line.pre.link_forward(current_new_line)
                    else:
                        start_base_line = current_new_line
                        start_base_e = current_height * p - a * (d - cross_d)
            e_diff = e_diff_last

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                current_new_line.link_forward(LinkedLinePDE(d, a))
                current_new_line = current_new_line.next

                # 上下半面 并没有赋值为1
                result[x_distance_index][y_distance_index] = 1
            else:
                result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从东向西过程
            x_grid_index -= 1
            x_distance_index -= 1

        # 最外层循环 对应网格点横坐标加1：从南向北过程
        # 拓展 x 边界
        start_distance_x_index += 1
        start_x_index += 1
        end_x_index -= 1

        y_grid_index += 1
        y_distance_index += 1
    # print("上面：\n", result)

    # 比较与原算法结果是否有异
    # true_result = judge_with_true.get_true_result()
    # print("是否相同：\n", judge_with_true.are_arrays_equal(result, true_result))

    u = 1 / (dem.dy * dem.rdy)

    # 左半面
    # 左半圆:从北往南，从东往西算
    is_center_right_top = dem.dy * lon_grid_center + dem.dx * lat_grid_center >= dem.dx * dem.dy

    start_index_adjust = 1 if is_center_right_top else 0
    end_index_adjust = -1 if is_center_right_bottom else 0

    start_y_index = y_grid_center + start_index_adjust + 1
    end_y_index = y_grid_center + end_index_adjust
    start_distance_y_index = y_grid_count + start_index_adjust

    x_grid_index = x_grid_center
    y_grid_index = start_y_index
    x_distance_index = x_grid_count - 1
    y_distance_index = start_distance_y_index
    last_height = dem.height[x_grid_index, y_grid_index + 1] - see_height
    start_base_line = None

    p = -1 / x_distance[x_distance_index]
    # 构造左半边初始参考线
    while y_grid_index >= end_y_index:
        # 当目标点的高度小于观察点的高度会存在负数的情况
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = -y_distance[y_distance_index] * p
        a = u * (current_height - last_height)

        if start_base_line is None:
            # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
            start_base_line = LinkedLinePDE(d, a)
            base_line = start_base_line
            start_base_e = current_height * p
        else:
            base_line.link_forward(LinkedLinePDE(d, a))
            base_line = base_line.next

        #  第一列是可见的
        result[x_distance_index][y_distance_index] = 1
        last_height = current_height

        y_grid_index -= 1
        y_distance_index -= 1

    start_distance_y_index += 1
    start_y_index += 1
    end_y_index -= 1
    x_grid_index -= 1
    x_distance_index -= 1

    # 左边逐层计算参考线与点的显隐性
    max_y_index = 2 * y_grid_count - 1
    while x_grid_index >= min_x:
        if start_distance_y_index > max_y_index:
            start_distance_y_index = max_y_index
            start_y_index -= 1
        y_grid_index = start_y_index
        y_distance_index = start_distance_y_index

        # 分析区域外拓展了一点
        last_height = dem.height[x_grid_index, y_grid_index + 1] - see_height

        p = -1 / x_distance[x_distance_index]
        base_line = start_base_line

        # 判断起点的显隐性
        # 初始化起点的d、a
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = -y_distance[y_distance_index] * p
        a = u * (current_height - last_height)

        while base_line.next.end_d < d:
            base_line = base_line.next
            # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
            start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        start_base_line = base_line
        current_base_e = start_base_e

        # 上一步相当于把延伸的点所参与的参考线段筛出去
        while base_line.end_d < d:
            base_line = base_line.next
            current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        """
        这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
        current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
        """
        e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

        # 可视性判断
        if e_diff >= 0:
            # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
            start_base_line = current_new_line = LinkedLinePDE(d, a)
            start_base_e = current_height * p
            result[x_distance_index][y_distance_index] = 1
        else:
            result[x_distance_index][y_distance_index] = 0
        last_d = d
        last_height = current_height

        # 从北向南过程
        y_grid_index -= 1
        y_distance_index -= 1

        # 判断后续点的显隐性别
        while y_grid_index >= min_y and y_grid_index >= end_y_index:
            # 形成新的调查线段
            d = -y_distance[y_distance_index] * p
            current_height = dem.height[x_grid_index, y_grid_index] - see_height
            a = u * (current_height - last_height)

            min_d = last_d
            # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
            while base_line.end_d < d:
                """
                求出轴斜率差值的增量
                这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                """
                e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                # 判断是否为交点 如果前后可视性不同则比存在交点
                # 当前是可视的
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = base_line.end_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last
                min_d = base_line.end_d
                base_line = base_line.next

            e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
            if e_diff >= 0:
                # 变为不可视
                if e_diff_last < 0:
                    # d很小时算不准
                    if e_diff < 5e-15:
                        cross_d = d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                    # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                    current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                    current_new_line = current_new_line.next

                    # 这条语句同时还更新了参考线
                    current_new_line.link_forward(base_line)
            else:
                if e_diff_last >= 0:
                    if e_diff > -5e-15:
                        cross_d = min_d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    current_new_line = LinkedLinePDE(cross_d, base_line.a)
                    if base_line.pre is not None:
                        base_line.pre.link_forward(current_new_line)
                    else:
                        start_base_line = current_new_line
                        start_base_e = current_height * p - a * (d - cross_d)
            e_diff = e_diff_last

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                current_new_line.link_forward(LinkedLinePDE(d, a))
                current_new_line = current_new_line.next
                result[x_distance_index][y_distance_index] = 1
            else:
                result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从北向南过程
            y_grid_index -= 1
            y_distance_index -= 1

        # 最外层循环 对应网格点横坐标加1：从东向西过程
        end_y_index -= 1
        start_distance_y_index += 1
        start_y_index += 1
        x_grid_index -= 1
        x_distance_index -= 1
    # print("左面：\n", result)

    # 上半面及下半面->横向格网间距实地距离的倒数 用它来计算线段系数a
    u = 1 / (dem.dx * dem.rdx)

    # 下半面
    # 下半圆：从北往南，从西往东算
    start_index_adjust = -1 if is_center_left_top else 0
    end_index_adjust = 1 if is_center_right_top else 0

    start_x_index = x_grid_center + start_index_adjust
    end_x_index = x_grid_center + end_index_adjust + 1
    start_distance_x_index = x_grid_count + start_index_adjust - 1

    x_grid_index = start_x_index
    y_grid_index = y_grid_center
    x_distance_index = start_distance_x_index
    y_distance_index = y_grid_count - 1
    last_height = dem.height[x_grid_index - 1, y_grid_index] - see_height
    start_base_line = None

    p = -1 / y_distance[y_distance_index]
    # 构造下半边初始参考线
    while x_grid_index <= end_x_index:
        # 当目标点的高度小于观察点的高度会存在负数的情况
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = x_distance[x_distance_index] * p
        a = u * (current_height - last_height)

        if start_base_line is None:
            # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
            start_base_line = LinkedLinePDE(d, a)
            base_line = start_base_line
            start_base_e = current_height * p
        else:
            base_line.link_forward(LinkedLinePDE(d, a))
            base_line = base_line.next

        #  第一列是可见的
        # result[x_distance_index][y_distance_index] = 1
        last_height = current_height

        x_grid_index += 1
        x_distance_index += 1

    start_x_index -= 1
    start_distance_x_index -= 1
    end_x_index += 1
    y_grid_index -= 1
    y_distance_index -= 1

    # 下边逐层计算参考线与点的显隐性
    while y_grid_index >= min_y:
        if start_distance_x_index < 0:
            start_distance_x_index = 0
            start_x_index += 1
        x_grid_index = start_x_index
        x_distance_index = start_distance_x_index

        # 分析区域外拓展了一点
        last_height = dem.height[x_grid_index - 1, y_grid_index] - see_height

        p = -1 / y_distance[y_distance_index]
        base_line = start_base_line

        # 判断起点的显隐性
        # 初始化起点的d、a
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = x_distance[x_distance_index] * p
        a = u * (current_height - last_height)

        """
        为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
        为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
        如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
        k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
        """
        while base_line.next.end_d < d:
            base_line = base_line.next
            # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
            start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        start_base_line = base_line
        current_base_e = start_base_e

        # 上一步相当于把延伸的点所参与的参考线段筛出去
        while base_line.end_d < d:
            base_line = base_line.next
            current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        """
        这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
        current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
        """
        e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

        # 可视性判断
        if e_diff >= 0:
            # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
            start_base_line = current_new_line = LinkedLinePDE(d, a)
            start_base_e = current_height * p
            # result[x_distance_index][y_distance_index] = 1
        else:
            result[x_distance_index][y_distance_index] = 0

        last_d = d
        last_height = current_height

        # 从西向东过程
        x_grid_index += 1
        x_distance_index += 1

        # 判断后续点的显隐性别
        while x_grid_index <= max_x and x_grid_index <= end_x_index:
            # 形成新的调查线段
            d = x_distance[x_distance_index] * p
            current_height = dem.height[x_grid_index, y_grid_index] - see_height
            a = u * (current_height - last_height)

            min_d = last_d
            # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
            while base_line.end_d < d:
                """
                求出轴斜率差值的增量
                这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                """
                e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                # 判断是否为交点 如果前后可视性不同则比存在交点
                # 当前是可视的
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = base_line.end_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last
                min_d = base_line.end_d
                base_line = base_line.next

            e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
            if e_diff >= 0:
                # 变为不可视
                if e_diff_last < 0:
                    # d很小时算不准
                    if e_diff < 5e-15:
                        cross_d = d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                    # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                    current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                    current_new_line = current_new_line.next

                    # 这条语句同时还更新了参考线
                    current_new_line.link_forward(base_line)
            else:
                if e_diff_last >= 0:
                    if e_diff > -5e-15:
                        cross_d = min_d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    current_new_line = LinkedLinePDE(cross_d, base_line.a)
                    if base_line.pre is not None:
                        base_line.pre.link_forward(current_new_line)
                    else:
                        start_base_line = current_new_line
                        start_base_e = current_height * p - a * (d - cross_d)
            e_diff = e_diff_last

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                current_new_line.link_forward(LinkedLinePDE(d, a))
                current_new_line = current_new_line.next

                result[x_distance_index][y_distance_index] = 1
            else:
                result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从西向东过程
            x_grid_index += 1
            x_distance_index += 1

        # 最外层循环 对应网格点横坐标加1：从北向南过程
        start_x_index -= 1
        start_distance_x_index -= 1
        end_x_index += 1
        y_grid_index -= 1
        y_distance_index -= 1

    # print("最终结果:\n", result)

    # print("最终是否相同：\n", judge_with_true.are_arrays_equal(result, true_result))
    return result


def analysis_by_xpderl_with_angle(x_center, y_center, h_center, radius, h_stand, horizontal_start_angle,
                                  horizontal_end_angle, current_dem=dem):
    # 初始化
    global dem
    dem = current_dem
    # 初始化
    x_grid_observe, y_grid_observe, dem_minx, dem_miny, x_grid_count, y_grid_count, see_height, x_grid_center, y_grid_center, min_x, max_x, min_y, max_y, is_in_grid_point, result = get_initial_param_spderl(
        x_center, y_center, h_center, radius, h_stand)
    x_distance, y_distance = get_xy_distance(x_center, y_center, dem_minx, dem_miny, x_grid_count, y_grid_count)

    # 纵向格网间距实地距离的倒数 用它来计算线段系数a
    u = 1 / (dem.dy * dem.rdy)

    # p,d根据不同区域会有所变化，但均符合以下两点：①p始终保持正值 ②d值随着调查线的构建应该逐渐变大
    # 邻近值p->1/x
    p = 0.0

    # 当前经纬度对应的高程值
    current_height = 0.0
    last_height = 0.0

    # 方向值d->y/x
    d = 0.0
    # 用于对比旧基准线上一次的当前点
    last_d = 0.0
    # 交点方向值
    cross_d = 0.0
    # 记录一次求交运算中交点区间的最小值
    min_d = 0.0

    # 线段对应系数a->u*(current_height-last_height)
    a = 0.0

    # 第一个点对应的e->height*p
    start_base_e = 0.0

    current_base_e = 0.0

    # 参考线的起始线段
    start_base_line = None
    # 用于对比旧基准线的当前点
    base_line = None
    # 即将接入新参考线的一部分线段,在不确定该线段定义域范围时暂存于该变量
    current_new_line = None

    # 求出中心点在中心格网的相对位置 这里的中心网格不代表为整个结果的中心点 而是代表对一个网格框(由四个点组成)他们的中心
    # dem_min对应为一个网格点的经纬度

    lon_grid_center = (x_center - dem_minx) % dem.dx
    lat_grid_center = (y_center - dem_miny) % dem.dy

    # 针对不同区域来判断
    # 观察点处于对角线上==>分割线与 DEM 网格中的对角线重合，则无需扩展计算区域
    # 如果与某一条对角线重合也不用扩展
    # 拓展边界应进行调整 若在对角线上则无需拓展 反之拓展
    is_center_left_top = lat_grid_center / lon_grid_center - dem.dy / dem.dx >= 0
    # 结合图像的几何意义更易于理解
    is_center_left_bottom = dem.dy * lon_grid_center + dem.dx * lat_grid_center <= dem.dx * dem.dy

    # 格网点序号值,用来获取对应坐标的高程值height=dem[x_grid_index,y_grid_index]
    x_grid_index = 0
    y_grid_index = 0

    # x、y轴距离索引,用来获取对应地心坐标系x、y值 x=x_distance[x_distance_index]->将其转换为p、d、e
    x_distance_index = 0
    y_distance_index = 0

    # 边界索引值
    x_right_index = 2 * x_grid_count - 1
    y_top_index = 2 * y_grid_count - 1

    # 当前点对应参考线和调查线的e差值 dif_e=current_e(调查线)-current_base_e(参考线)
    e_diff = 0.0
    e_diff_last = 0.0

    # 右半面
    # 以下算法认为横向间隔和纵向间隔相同
    # 右半圆: 从南往北,从西往东算
    # 通过约束纵向y索引来限制判断区域
    # 拓展边界
    start_index_adjust = -1 if is_center_left_bottom else 0
    end_index_adjust = 1 if is_center_left_top else 0

    # 约束循环条件
    start_y_index = y_grid_center + start_index_adjust
    end_y_index = y_grid_center + end_index_adjust + 1

    # 调整对应的d边界
    start_distance_y_index = y_grid_count + start_index_adjust - 1

    x_grid_index = x_grid_center + 1
    y_grid_index = start_y_index
    x_distance_index = x_grid_count
    y_distance_index = start_distance_y_index
    last_height = dem.height[x_grid_index, y_grid_index - 1] - see_height
    p = 1 / x_distance[x_distance_index]

    # 观察点不是网格点该算法也只是将其视为中心网格点上方或者下方的网格点 并没有对其求值
    # 构造右半边初始参考线
    while y_grid_index <= end_y_index:
        # 当目标点的高度小于观察点的高度会存在负数的情况
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = y_distance[y_distance_index] * p
        a = u * (current_height - last_height)

        if start_base_line is None:
            # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
            start_base_line = LinkedLinePDE(d, a)
            base_line = start_base_line
            start_base_e = current_height * p
        else:
            base_line.link_forward(LinkedLinePDE(d, a))
            base_line = base_line.next

        #  第一列是可见的
        result[x_distance_index][y_distance_index] = 1
        last_height = current_height

        y_grid_index += 1
        y_distance_index += 1

    # 区域内第一个观察点的y坐标
    start_distance_y_index -= 1

    # 向y方向上下拓展边界 x型分区
    start_y_index -= 1
    end_y_index += 1

    # 目标点继续向右移动
    x_grid_index += 1
    x_distance_index += 1

    # 右边逐层计算参考线与点的显隐性
    while x_grid_index <= max_x:
        # 防止超出结果数组边界,回调索引
        if start_distance_y_index < 0:
            start_y_index += 1
            start_distance_y_index = 0

        y_grid_index = start_y_index
        y_distance_index = start_distance_y_index

        # 分析区域外拓展了一点
        last_height = dem.height[x_grid_index, y_grid_index - 1] - see_height

        p = 1 / x_distance[x_distance_index]
        base_line = start_base_line

        # 判断起点的显隐性
        # 初始化起点的d、a
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = y_distance[y_distance_index] * p
        a = u * (current_height - last_height)

        """
        为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
        为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
        如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
        k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
        """
        while base_line.next.end_d < d:
            base_line = base_line.next
            # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
            start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        start_base_line = base_line
        current_base_e = start_base_e

        # 上一步相当于把延伸的点所参与的参考线段筛出去
        while base_line.end_d < d:
            base_line = base_line.next
            current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        """
        这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
        current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
        """
        e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

        # 可视性判断
        if e_diff >= 0:
            # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
            start_base_line = current_new_line = LinkedLinePDE(d, a)
            start_base_e = current_height * p
            current_quadrant = 1 if y_grid_index >= y_grid_center else 4
            is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                            y_grid_observe, x_grid_index, y_grid_index, current_quadrant)
            if is_in_range and x_distance_index >= x_grid_count:
                result[x_distance_index][y_distance_index] = 1

        else:
            result[x_distance_index][y_distance_index] = 0
        last_d = d
        last_height = current_height

        # 从南向北过程
        y_grid_index += 1
        y_distance_index += 1

        # 判断后续点的显隐性别
        while y_grid_index <= max_y and y_grid_index <= end_y_index:
            # 形成新的调查线段
            d = y_distance[y_distance_index] * p
            current_height = dem.height[x_grid_index, y_grid_index] - see_height
            a = u * (current_height - last_height)

            min_d = last_d
            # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
            while base_line.end_d < d:
                """
                求出轴斜率差值的增量
                这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                """
                e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                # 判断是否为交点 如果前后可视性不同则比存在交点
                # 当前是可视的
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = base_line.end_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last
                min_d = base_line.end_d
                if base_line.next is None:
                    break
                base_line = base_line.next

            e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
            if e_diff >= 0:
                # 变为不可视
                if e_diff_last < 0:
                    # d很小时算不准
                    if e_diff < 5e-15:
                        cross_d = d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                    # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                    current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                    current_new_line = current_new_line.next

                    # 这条语句同时还更新了参考线
                    current_new_line.link_forward(base_line)
            else:
                if e_diff_last >= 0:
                    if e_diff > -5e-15:
                        cross_d = min_d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    current_new_line = LinkedLinePDE(cross_d, base_line.a)
                    if base_line.pre is not None:
                        base_line.pre.link_forward(current_new_line)
                    else:
                        start_base_line = current_new_line
                        start_base_e = current_height * p - a * (d - cross_d)
            e_diff = e_diff_last

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                current_new_line.link_forward(LinkedLinePDE(d, a))
                current_new_line = current_new_line.next
                result[x_distance_index][y_distance_index] = 1
            else:
                result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从南向北过程
            y_grid_index += 1
            y_distance_index += 1

        # 最外层循环 对应网格点横坐标加1：从西向东过程
        x_grid_index += 1
        x_distance_index += 1
        start_y_index -= 1
        start_distance_y_index -= 1
        end_y_index += 1
    # print("右面：\n", result)

    # 上半面及下半面->横向格网间距实地距离的倒数 用它来计算线段系数a
    u = 1 / (dem.dx * dem.rdx)

    # 上半面
    # 上半圆:从南往北，从东往西算
    is_center_right_bottom = lat_grid_center / lon_grid_center - dem.dy / dem.dx <= 0

    start_index_adjust = 1 if is_center_right_bottom else 0
    end_index_adjust = -1 if is_center_left_bottom else 0

    start_x_index = x_grid_center + start_index_adjust + 1
    start_distance_x_index = x_grid_count + start_index_adjust
    end_x_index = x_grid_center + end_index_adjust

    x_grid_index = start_x_index
    y_grid_index = y_grid_center + 1
    x_distance_index = start_distance_x_index
    y_distance_index = y_grid_count
    last_height = dem.height[x_grid_index + 1, y_grid_index] - see_height
    start_base_line = None

    p = 1 / y_distance[y_distance_index]

    # 构造上半边初始参考线
    while x_grid_index >= end_x_index:
        # 当目标点的高度小于观察点的高度会存在负数的情况
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = -x_distance[x_distance_index] * p
        a = u * (current_height - last_height)

        if start_base_line is None:
            # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
            start_base_line = LinkedLinePDE(d, a)
            base_line = start_base_line
            start_base_e = current_height * p
        else:
            base_line.link_forward(LinkedLinePDE(d, a))
            base_line = base_line.next

        #  第一列是可见的
        # result[x_distance_index][y_distance_index] = 1
        last_height = current_height

        x_grid_index -= 1
        x_distance_index -= 1

    # 拓展 x 方向
    start_x_index += 1
    start_distance_x_index += 1
    end_x_index -= 1

    y_grid_index += 1
    y_distance_index += 1

    # 上边逐层计算参考线与点的显隐性
    # 防止超出x边界
    max_right_index = 2 * x_grid_count - 1
    while y_grid_index <= max_y:
        if start_distance_x_index > max_right_index:
            start_x_index -= 1
            start_distance_x_index = max_right_index
        x_grid_index = start_x_index
        x_distance_index = start_distance_x_index

        # 分析区域外拓展了一点
        last_height = dem.height[x_grid_index + 1, y_grid_index] - see_height

        p = 1 / y_distance[y_distance_index]
        base_line = start_base_line

        # 判断起点的显隐性
        # 初始化起点的d、a
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = -x_distance[x_distance_index] * p
        a = u * (current_height - last_height)

        """
        为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
        为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
        如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
        k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
        """
        while base_line.next.end_d < d:
            base_line = base_line.next
            # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
            start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        start_base_line = base_line
        current_base_e = start_base_e

        # 上一步相当于把延伸的点所参与的参考线段筛出去
        while base_line.end_d < d:
            base_line = base_line.next
            current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        """
        这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
        current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
        """
        e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

        # 可视性判断
        if e_diff >= 0:
            # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
            start_base_line = current_new_line = LinkedLinePDE(d, a)
            start_base_e = current_height * p
            # result[x_distance_index][y_distance_index] = 1
        else:
            result[x_distance_index][y_distance_index] = 0
        last_d = d
        last_height = current_height

        # 从东向西过程
        x_grid_index -= 1
        x_distance_index -= 1

        # 判断后续点的显隐性别
        while x_grid_index >= min_x and x_grid_index >= end_x_index:
            # 形成新的调查线段
            d = -x_distance[x_distance_index] * p
            current_height = dem.height[x_grid_index, y_grid_index] - see_height
            a = u * (current_height - last_height)

            min_d = last_d
            # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
            while base_line.end_d < d:
                """
                求出轴斜率差值的增量
                这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                """
                e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                # 判断是否为交点 如果前后可视性不同则比存在交点
                # 当前是可视的
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = base_line.end_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last
                min_d = base_line.end_d
                if base_line.next is None:
                    break
                base_line = base_line.next

            e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
            if e_diff >= 0:
                # 变为不可视
                if e_diff_last < 0:
                    # d很小时算不准
                    if e_diff < 5e-15:
                        cross_d = d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                    # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                    current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                    current_new_line = current_new_line.next

                    # 这条语句同时还更新了参考线
                    current_new_line.link_forward(base_line)
            else:
                if e_diff_last >= 0:
                    if e_diff > -5e-15:
                        cross_d = min_d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    current_new_line = LinkedLinePDE(cross_d, base_line.a)
                    if base_line.pre is not None:
                        base_line.pre.link_forward(current_new_line)
                    else:
                        start_base_line = current_new_line
                        start_base_e = current_height * p - a * (d - cross_d)
            e_diff = e_diff_last

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                current_new_line.link_forward(LinkedLinePDE(d, a))
                current_new_line = current_new_line.next
                current_quadrant = 1 if x_grid_index >= x_grid_center else 2
                is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                                y_grid_observe, x_grid_index, y_grid_index, current_quadrant)
                if is_in_range and y_distance_index >= y_grid_count:
                    result[x_distance_index][y_distance_index] = 1
            else:
                result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从东向西过程
            x_grid_index -= 1
            x_distance_index -= 1

        # 最外层循环 对应网格点横坐标加1：从南向北过程
        # 拓展 x 边界
        start_distance_x_index += 1
        start_x_index += 1
        end_x_index -= 1

        y_grid_index += 1
        y_distance_index += 1
    # print("上面：\n", result)

    u = 1 / (dem.dy * dem.rdy)

    # 左半面
    # 左半圆:从北往南，从东往西算
    is_center_right_top = dem.dy * lon_grid_center + dem.dx * lat_grid_center >= dem.dx * dem.dy

    start_index_adjust = 1 if is_center_right_top else 0
    end_index_adjust = -1 if is_center_right_bottom else 0

    start_y_index = y_grid_center + start_index_adjust + 1
    end_y_index = y_grid_center + end_index_adjust
    start_distance_y_index = y_grid_count + start_index_adjust

    x_grid_index = x_grid_center
    y_grid_index = start_y_index
    x_distance_index = x_grid_count - 1
    y_distance_index = start_distance_y_index
    last_height = dem.height[x_grid_index, y_grid_index + 1] - see_height
    start_base_line = None

    p = -1 / x_distance[x_distance_index]
    # 构造左半边初始参考线
    while y_grid_index >= end_y_index:
        # 当目标点的高度小于观察点的高度会存在负数的情况
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = -y_distance[y_distance_index] * p
        a = u * (current_height - last_height)

        if start_base_line is None:
            # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
            start_base_line = LinkedLinePDE(d, a)
            base_line = start_base_line
            start_base_e = current_height * p
        else:
            base_line.link_forward(LinkedLinePDE(d, a))
            base_line = base_line.next

        #  第一列是可见的
        result[x_distance_index][y_distance_index] = 1
        last_height = current_height

        y_grid_index -= 1
        y_distance_index -= 1

    start_distance_y_index += 1
    start_y_index += 1
    end_y_index -= 1
    x_grid_index -= 1
    x_distance_index -= 1

    # 左边逐层计算参考线与点的显隐性
    max_y_index = 2 * y_grid_count - 1
    while x_grid_index >= min_x:
        if start_distance_y_index > max_y_index:
            start_distance_y_index = max_y_index
            start_y_index -= 1
        y_grid_index = start_y_index
        y_distance_index = start_distance_y_index

        # 分析区域外拓展了一点
        last_height = dem.height[x_grid_index, y_grid_index + 1] - see_height

        p = -1 / x_distance[x_distance_index]
        base_line = start_base_line

        # 判断起点的显隐性
        # 初始化起点的d、a
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = -y_distance[y_distance_index] * p
        a = u * (current_height - last_height)

        while base_line.next.end_d < d:
            base_line = base_line.next
            # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
            start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        start_base_line = base_line
        current_base_e = start_base_e

        # 上一步相当于把延伸的点所参与的参考线段筛出去
        while base_line.end_d < d:
            base_line = base_line.next
            current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        """
        这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
        current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
        """
        e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

        # 可视性判断
        if e_diff >= 0:
            # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
            start_base_line = current_new_line = LinkedLinePDE(d, a)
            start_base_e = current_height * p
            current_quadrant = 2 if y_grid_index >= y_grid_center else 3
            is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                            y_grid_observe, x_grid_index, y_grid_index, current_quadrant)
            if is_in_range and x_distance_index <= x_grid_count:
                result[x_distance_index][y_distance_index] = 1
        else:
            result[x_distance_index][y_distance_index] = 0
        last_d = d
        last_height = current_height

        # 从北向南过程
        y_grid_index -= 1
        y_distance_index -= 1

        # 判断后续点的显隐性别
        while y_grid_index >= min_y and y_grid_index >= end_y_index:
            # 形成新的调查线段
            d = -y_distance[y_distance_index] * p
            current_height = dem.height[x_grid_index, y_grid_index] - see_height
            a = u * (current_height - last_height)

            min_d = last_d
            # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
            while base_line.end_d < d:
                """
                求出轴斜率差值的增量
                这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                """
                e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                # 判断是否为交点 如果前后可视性不同则比存在交点
                # 当前是可视的
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = base_line.end_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last
                min_d = base_line.end_d
                if base_line.next is None:
                    break
                base_line = base_line.next

            e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
            if e_diff >= 0:
                # 变为不可视
                if e_diff_last < 0:
                    # d很小时算不准
                    if e_diff < 5e-15:
                        cross_d = d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                    # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                    current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                    current_new_line = current_new_line.next

                    # 这条语句同时还更新了参考线
                    current_new_line.link_forward(base_line)
            else:
                if e_diff_last >= 0:
                    if e_diff > -5e-15:
                        cross_d = min_d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    current_new_line = LinkedLinePDE(cross_d, base_line.a)
                    if base_line.pre is not None:
                        base_line.pre.link_forward(current_new_line)
                    else:
                        start_base_line = current_new_line
                        start_base_e = current_height * p - a * (d - cross_d)
            e_diff = e_diff_last

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                current_new_line.link_forward(LinkedLinePDE(d, a))
                current_new_line = current_new_line.next
                current_quadrant = 2 if y_grid_index >= y_grid_center else 3
                is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                                y_grid_observe, x_grid_index, y_grid_index, current_quadrant)
                if is_in_range and x_distance_index <= x_grid_count:
                    result[x_distance_index][y_distance_index] = 1
            else:
                result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从北向南过程
            y_grid_index -= 1
            y_distance_index -= 1

        # 最外层循环 对应网格点横坐标加1：从东向西过程
        end_y_index -= 1
        start_distance_y_index += 1
        start_y_index += 1
        x_grid_index -= 1
        x_distance_index -= 1
    # print("左面：\n", result)

    # 上半面及下半面->横向格网间距实地距离的倒数 用它来计算线段系数a
    u = 1 / (dem.dx * dem.rdx)

    # 下半面
    # 下半圆：从北往南，从西往东算
    start_index_adjust = -1 if is_center_left_top else 0
    end_index_adjust = 1 if is_center_right_top else 0

    start_x_index = x_grid_center + start_index_adjust
    end_x_index = x_grid_center + end_index_adjust + 1
    start_distance_x_index = x_grid_count + start_index_adjust - 1

    x_grid_index = start_x_index
    y_grid_index = y_grid_center
    x_distance_index = start_distance_x_index
    y_distance_index = y_grid_count - 1
    last_height = dem.height[x_grid_index - 1, y_grid_index] - see_height
    start_base_line = None

    p = -1 / y_distance[y_distance_index]
    # 构造下半边初始参考线
    while x_grid_index <= end_x_index:
        # 当目标点的高度小于观察点的高度会存在负数的情况
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = x_distance[x_distance_index] * p
        a = u * (current_height - last_height)

        if start_base_line is None:
            # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
            start_base_line = LinkedLinePDE(d, a)
            base_line = start_base_line
            start_base_e = current_height * p
        else:
            base_line.link_forward(LinkedLinePDE(d, a))
            base_line = base_line.next

        #  第一列是可见的
        # result[x_distance_index][y_distance_index] = 1
        last_height = current_height

        x_grid_index += 1
        x_distance_index += 1

    start_x_index -= 1
    start_distance_x_index -= 1
    end_x_index += 1
    y_grid_index -= 1
    y_distance_index -= 1

    # 下边逐层计算参考线与点的显隐性
    while y_grid_index >= min_y:
        if start_distance_x_index < 0:
            start_distance_x_index = 0
            start_x_index += 1
        x_grid_index = start_x_index
        x_distance_index = start_distance_x_index

        # 分析区域外拓展了一点
        last_height = dem.height[x_grid_index - 1, y_grid_index] - see_height

        p = -1 / y_distance[y_distance_index]
        base_line = start_base_line

        # 判断起点的显隐性
        # 初始化起点的d、a
        current_height = dem.height[x_grid_index, y_grid_index] - see_height
        d = x_distance[x_distance_index] * p
        a = u * (current_height - last_height)

        """
        为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
        为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
        如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
        k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
        """
        while base_line.next.end_d < d:
            base_line = base_line.next
            # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
            start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        start_base_line = base_line
        current_base_e = start_base_e

        # 上一步相当于把延伸的点所参与的参考线段筛出去
        while base_line.end_d < d:
            base_line = base_line.next
            current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

        """
        这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
        current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
        """
        e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

        # 可视性判断
        if e_diff >= 0:
            # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
            start_base_line = current_new_line = LinkedLinePDE(d, a)
            start_base_e = current_height * p
            # result[x_distance_index][y_distance_index] = 1
        else:
            result[x_distance_index][y_distance_index] = 0

        last_d = d
        last_height = current_height

        # 从西向东过程
        x_grid_index += 1
        x_distance_index += 1

        # 判断后续点的显隐性别
        while x_grid_index <= max_x and x_grid_index <= end_x_index:
            # 形成新的调查线段
            d = x_distance[x_distance_index] * p
            current_height = dem.height[x_grid_index, y_grid_index] - see_height
            a = u * (current_height - last_height)

            min_d = last_d
            # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
            while base_line.end_d < d:
                """
                求出轴斜率差值的增量
                这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                """
                e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                # 判断是否为交点 如果前后可视性不同则比存在交点
                # 当前是可视的
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = base_line.end_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last
                min_d = base_line.end_d
                if base_line.next is None:
                    break
                base_line = base_line.next

            e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
            if e_diff >= 0:
                # 变为不可视
                if e_diff_last < 0:
                    # d很小时算不准
                    if e_diff < 5e-15:
                        cross_d = d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                    # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                    current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                    current_new_line = current_new_line.next

                    # 这条语句同时还更新了参考线
                    current_new_line.link_forward(base_line)
            else:
                if e_diff_last >= 0:
                    if e_diff > -5e-15:
                        cross_d = min_d
                    else:
                        # 交点的d值
                        cross_d = min_d + e_diff / (base_line.a - a)

                    current_new_line = LinkedLinePDE(cross_d, base_line.a)
                    if base_line.pre is not None:
                        base_line.pre.link_forward(current_new_line)
                    else:
                        start_base_line = current_new_line
                        start_base_e = current_height * p - a * (d - cross_d)
            e_diff = e_diff_last

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                current_new_line.link_forward(LinkedLinePDE(d, a))
                current_new_line = current_new_line.next
                current_quadrant = 4 if x_grid_index >= x_grid_center else 3
                is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                                y_grid_observe, x_grid_index, y_grid_index, current_quadrant)
                if is_in_range and y_distance_index <= y_grid_count:
                    result[x_distance_index][y_distance_index] = 1
            else:
                result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从西向东过程
            x_grid_index += 1
            x_distance_index += 1

        # 最外层循环 对应网格点横坐标加1：从北向南过程
        start_x_index -= 1
        start_distance_x_index -= 1
        end_x_index += 1
        y_grid_index -= 1
        y_distance_index -= 1

    # print("最终结果:\n", result)
    return result


# 扇形视域分析
"""
r: 视距
horizontal_start_angle: 水平扫描起始角度0-360 其中0指向北
horizontal_end_angle： 水平扫描截至角度0-360
"""


def analysis_by_dpderl_simplified(x_center, y_center, h_center, radius, h_stand, horizontal_start_angle=0.0,
                                  horizontal_end_angle=360.0, current_dem=dem):
    global dem
    dem = current_dem
    # 初始化
    # 其中x_grid_count --> x_grad_center + 1、 min_y --> 0、 max_y --> 2*y_grid_count - 1
    # 对于dem.height 属于[min_x,max_x],[min_y,max_y] result属于[0,2*x_grid_count-1]
    x_grid_observe, y_grid_observe, dem_minx, dem_miny, x_grid_count, y_grid_count, see_height, x_grid_center, y_grid_center, min_x, max_x, min_y, max_y, is_in_grid_point, result = get_initial_param_spderl(
        x_center, y_center, h_center, radius, h_stand)
    x_distance, y_distance = get_xy_distance(x_center, y_center, dem_minx, dem_miny, x_grid_count, y_grid_count)

    result_spderl = np.zeros((2 * x_grid_count, 2 * y_grid_count))
    result_opt_rl_spderl = np.zeros((2 * x_grid_count, 2 * y_grid_count))
    print(
        f"x_center:{x_center},y_center:{y_center}\nx_grid_observe:{x_grid_observe}\ny_grid_observe:{y_grid_observe}\nx_grid_count:{x_grid_count}\ny_grid_count:{y_grid_count}\nx_grid_center:{x_grid_center}\ny_grid_center:{y_grid_center}\nmin_x:{min_x}\nmax_x:{max_x}\nmin_y:{min_y}\nmax_y:{max_y}\n")
    # 实际区域格网点数
    count = 0

    # 纵向格网间距实地距离的倒数 用它来计算线段系数a
    u = 1 / (dem.dy * dem.rdy)

    # 第一个点对应的e->height*p
    start_base_e = 0.0
    # 参考线的起始线段
    start_base_line = None
    # 用于对比旧基准线的当前点
    base_line = None
    # 即将接入新参考线的一部分线段,在不确定该线段定义域范围时暂存于该变量
    current_new_line = None

    # 观察者对面对的网格线方向必须为横线才可以
    # 此处角度分析不考虑最终角度大于360°的情况 且终止角度>起始角度
    # 根据角度值来划分区域并判断各区域的方向
    angle_range = horizontal_end_angle - horizontal_start_angle

    # 求出中心点在中心格网的相对位置 这里的中心网格不代表为整个结果的中心点 而是代表对一个网格框(由四个点组成)他们的中心
    # dem_min对应为一个网格点的经纬度
    lon_grid_center = (x_center - dem_minx) % dem.dx
    lat_grid_center = (y_center - dem_miny) % dem.dy

    # 判断观察点是否位于左上或者右下对角线,根据具体分析区域判断加还是减
    is_observe_not_in_left_top_diagonal = 1 if lat_grid_center / lon_grid_center != dem.dy / dem.dx else 0
    is_observe_not_in_right_bottom_diagonal = 1 if dem.dy * lon_grid_center + dem.dx * lat_grid_center != dem.dx * dem.dy else 0
    partition_strategy = PartitionAlgorithm(horizontal_start_angle, horizontal_end_angle,
                                            x_grid_observe, y_grid_observe, x_grid_count, y_grid_count, see_height,
                                            x_grid_center,
                                            y_grid_center, min_x, max_x, min_y, max_y, x_distance, y_distance, dem,
                                            is_observe_not_in_left_top_diagonal,
                                            is_observe_not_in_right_bottom_diagonal, result_spderl, quadrant,
                                            is_in_grid_point)
    partition_optimal_rl = PartitionOptimalRLAlgorithm(horizontal_start_angle, horizontal_end_angle,
                                                       x_grid_observe, y_grid_observe, x_grid_count, y_grid_count,
                                                       see_height,
                                                       x_grid_center,
                                                       y_grid_center, min_x, max_x, min_y, max_y, x_distance,
                                                       y_distance, dem,
                                                       is_observe_not_in_left_top_diagonal,
                                                       is_observe_not_in_right_bottom_diagonal, result_opt_rl_spderl,
                                                       quadrant,
                                                       is_in_grid_point)
    if angle_range <= 90:
        # 分析区域位于右侧 --> 纵向网格线(从南到北 从西到东)
        if horizontal_start_angle < 90:
            result_spderl, count = partition_strategy.analysis_right()
            result_opt_rl_spderl, _ = partition_optimal_rl.analysis_right()
        # 位于下侧 --> 横向网格线(从北向南，从西到东)
        elif 90 <= horizontal_start_angle < 180:
            result_spderl, count = partition_strategy.analysis_bottom()
            result_opt_rl_spderl, _ = partition_optimal_rl.analysis_bottom()

        # 位于左侧 --> 纵向网格线(从北向南,从东向西)
        elif 180 <= horizontal_start_angle < 270:
            result_spderl, count = partition_strategy.analysis_left()
            result_opt_rl_spderl, _ = partition_optimal_rl.analysis_left()

        # 位于上侧 --> 横向网格线(从南到北，从东到西)
        elif 270 <= horizontal_start_angle:
            result_spderl, count = partition_strategy.analysis_top()
            result_opt_rl_spderl, _ = partition_optimal_rl.analysis_top()
    elif 90 < angle_range < 180:
        # 存在两种分法
        # 分为两类：①大于180部分位于观察点下侧-->横向网格线
        #         ②小于180部分位于观察点右侧-->纵向网格线
        # 注意调中间索引
        # 通过三条直线方程来分割区域为下侧区域和右侧区域：l1->左下、l2->中间区域、l3->右边
        if horizontal_start_angle < 90 and 180 <= horizontal_end_angle:
            result_spderl, count = partition_strategy.analysis_right_bottom()
            result_opt_rl_spderl, _ = partition_optimal_rl.analysis_right_bottom()

        # 分为两类：①小于180部分位于观察点下侧-->横向网格线
        #         ②大于180部分位于观察点左侧-->纵向网格线
        # 通过三条直线方程来分割区域为下侧区域和左侧区域：l1->左上、l2->中间区域、l3->右下
        elif 90 < horizontal_start_angle < 180 and 270 <= horizontal_end_angle:
            result_spderl, count = partition_strategy.analysis_left_bottom()
            result_opt_rl_spderl, _ = partition_optimal_rl.analysis_left_bottom()

        # 位于右侧
        elif horizontal_start_angle < 90 and horizontal_end_angle < 180:
            result_spderl, count = partition_strategy.analysis_right()
            result_opt_rl_spderl, _ = partition_optimal_rl.analysis_right()

        # 位于左侧
        elif 180 <= horizontal_start_angle and horizontal_end_angle < 360:
            result_spderl, count = partition_strategy.analysis_left()
            result_opt_rl_spderl, _ = partition_optimal_rl.analysis_left()

        # 位于下侧
        elif 90 < horizontal_start_angle and horizontal_end_angle < 270:
            result_spderl, count = partition_strategy.analysis_bottom()
            result_opt_rl_spderl, _ = partition_optimal_rl.analysis_bottom()
    elif 180 <= angle_range < 270:
        # 分为两类：①大于180部分位于观察点下侧-->横向网格线
        #         ②小于180部分位于观察点右侧-->纵向网格线
        # 若终止角为180则表面没有下侧 只有右侧
        # 通过三条直线方程来分割区域为下侧区域和右侧区域：l1->左下、l2->中间区域、l3->右边
        if horizontal_start_angle < 90 and 180 <= horizontal_end_angle < 270:
            result_spderl, count = partition_strategy.analysis_right_bottom()
            result_opt_rl_spderl, _ = partition_optimal_rl.analysis_right_bottom()

        # 分为两类：①大于180部分位于观察点左侧-->纵向网格线
        #         ②小于180部分位于观察点右侧-->纵向网格线
        # 通过三条直线方程来分割区域为左侧区域和右侧区域：l1->左、l2->中间区域、l3->右
        elif horizontal_start_angle < 90 and 270 <= horizontal_end_angle:
            result_spderl, count = partition_strategy.analysis_left_right()
            result_opt_rl_spderl, _ = partition_optimal_rl.analysis_left_right()

        # 分为两类：①大于180部分位于观察点左侧-->纵向网格线
        #         ②小于180部分位于观察点右侧-->纵向网格线 ×
        #         ③若大于90度分析区域位于观察点下侧此时因采用横向网格线 √
        # 通过三条直线方程来分割区域为左侧区域和右侧区域：l1->左、l2->中间区域、l3->右
        elif 90 <= horizontal_start_angle and 270 <= horizontal_end_angle:
            result_spderl, count = partition_strategy.analysis_left_bottom()
            result_opt_rl_spderl, _ = partition_optimal_rl.analysis_left_bottom()
    elif 270 <= angle_range < 360:
        result_spderl, count = partition_strategy.analysis_left_right()
        result_opt_rl_spderl, _ = partition_optimal_rl.analysis_left_right()
    else:
        # 与XPDERL算法一致
        to_x = x_center + radius / dem.rdx
        to_y = y_center
        a = analysis_by_xpderl(x_center, y_center, h_center, to_x, to_y, h_stand)
        return a, a.size
    return result_spderl, count, result_opt_rl_spderl


def analysis_by_r3(x_center, y_center, h_center, radius, h_stand, horizontal_start_angle,
                   horizontal_end_angle, current_dem=dem):
    global dem
    dem = current_dem
    # 初始化
    x_grid_observe, y_grid_observe, dem_minx, dem_miny, x_grid_count, y_grid_count, see_height, x_grid_center, y_grid_center, min_x, max_x, min_y, max_y, is_in_grid_point, result = get_initial_param_spderl(
        x_center, y_center, h_center, radius, h_stand)
    x_distance, y_distance = get_xy_distance(x_center, y_center, dem_minx, dem_miny, x_grid_count, y_grid_count)

    width = 2 * x_grid_count
    height = 2 * y_grid_count
    real_lon = dem.real_distance
    real_lat = dem.real_distance
    # 当前所在象限
    current_quadrant = 1

    # 第一象限
    i = x_grid_count
    while i < width:
        j = y_grid_count
        while j < height:
            # 观察点到目标点的距离
            r = math.sqrt(x_distance[i] ** 2 + y_distance[j] ** 2)
            # 高度比距离，观察点至目标点仰角
            s = (dem.height[min_x + i][min_y + j] - see_height) / r
            max_s = s
            # 当前观察点至目标点与横向格网线的交点数
            x_intersection = math.floor(x_distance[i] / real_lon)

            # 获取当前目标点单个网格 LOS对应的dx[i] 的dy[j]
            dp_lat = y_distance[j] / x_distance[i] * real_lon

            # 格网点对角线长度
            dp_l = math.sqrt(dp_lat ** 2 + real_lon ** 2)

            # 循环遍历x交点，看纵向网格线，从目标点到观察点进行计算的

            k = 1
            while k <= x_intersection:
                # 第k个交点时对应的y长度
                dy = dp_lat * k

                # 交点至观察点线的长度
                length = r - dp_l * k

                lon_index = i - k + min_x

                # 表示目标点在当前纬度格点之间的相对位置
                inner_y = dy % real_lat
                # 向下取整 负数的时候会与int强转有差异 -2.1 --> -3 而不是-2
                lat_index = j - math.floor(dy / real_lat) + min_y

                # 线性插值
                h = inner_y * (dem.height[lon_index][lat_index - 1] - dem.height[lon_index][lat_index]) / real_lat + \
                    dem.height[lon_index][lat_index] - see_height

                current_s = h / length

                if current_s > max_s:
                    max_s = current_s
                k += 1

            y_intersection = math.floor(y_distance[j] / real_lat)
            dp_lon = x_distance[i] / y_distance[j] * real_lat
            dp_l = math.sqrt(dp_lon ** 2 + real_lat ** 2)
            k = 1
            while k <= y_intersection:
                dx = dp_lon * k
                length = r - dp_l * k

                lat_index = j - k + min_y
                inner_x = dx % real_lon
                lon_index = i - math.floor(dx / real_lon) + min_x

                h = inner_x * (dem.height[lon_index - 1][lat_index] - dem.height[lon_index][lat_index]) / real_lon + \
                    dem.height[lon_index][lat_index] - see_height

                current_s = h / length
                if current_s > max_s:
                    max_s = current_s
                k += 1

            is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                            y_grid_observe, min_x + i, min_y + j, current_quadrant)
            if s >= max_s and is_in_range:
                result[i][j] = 1
            else:
                result[i][j] = 0
            j += 1
        i += 1
    # (f"第一象限：\n{result}")

    current_quadrant = 2
    # 第二象限
    i = 0
    while i < x_grid_count:
        j = y_grid_count
        while j < height:
            # 观察点到目标点的距离
            r = math.sqrt(x_distance[i] ** 2 + y_distance[j] ** 2)
            # 高度比距离，观察点至目标点仰角
            s = (dem.height[min_x + i][min_y + j] - see_height) / r
            max_s = s

            # 当前观察点至目标点与横向格网线的交点数
            x_intersection = math.floor(-x_distance[i] / real_lon)

            # 获取当前目标点单个网格 LOS对应的dx[i] 的dy[j]
            dp_lat = y_distance[j] / -x_distance[i] * real_lon

            # 格网点对角线长度
            dp_l = math.sqrt(dp_lat ** 2 + real_lon ** 2)

            # 循环遍历x交点，看纵向网格线，从目标点到观察点进行计算的

            k = 1
            while k <= x_intersection:
                # 第k个交点时对应的y长度
                dy = dp_lat * k

                # 交点至观察点线的长度
                length = r - dp_l * k

                lon_index = i + k + min_x

                # 表示目标点在当前纬度格点之间的相对位置
                inner_y = dy % real_lat
                lat_index = j - math.floor(dy / real_lat) + min_y

                # 线性插值
                h = inner_y * (dem.height[lon_index][lat_index - 1] - dem.height[lon_index][lat_index]) / real_lat + \
                    dem.height[lon_index][lat_index] - see_height

                current_s = h / length

                if current_s > max_s:
                    max_s = current_s
                k += 1

            y_intersection = math.floor(y_distance[j] / real_lat)
            dp_lon = -x_distance[i] / y_distance[j] * real_lat
            dp_l = math.sqrt(dp_lon ** 2 + real_lat ** 2)
            k = 1
            while k <= y_intersection:
                dx = dp_lon * k
                length = r - dp_l * k

                lat_index = j - k + min_y
                inner_x = dx % real_lon
                lon_index = i + math.floor(dx / real_lon) + min_x

                h = inner_x * (dem.height[lon_index + 1][lat_index] - dem.height[lon_index][lat_index]) / real_lon + \
                    dem.height[lon_index][lat_index] - see_height

                current_s = h / length
                if current_s > max_s:
                    max_s = current_s
                k += 1

            is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                            y_grid_observe, min_x + i, min_y + j, current_quadrant)
            if s >= max_s and is_in_range:
                result[i][j] = 1
            else:
                result[i][j] = 0
            j += 1
        i += 1
    # print(f"第二象限：\n{result}")

    current_quadrant = 3
    # 第三象限
    i = 0
    while i < x_grid_count:
        j = 0
        while j < y_grid_count:
            # 观察点到目标点的距离
            r = math.sqrt(x_distance[i] ** 2 + y_distance[j] ** 2)
            # 高度比距离，观察点至目标点仰角
            s = (dem.height[min_x + i][min_y + j] - see_height) / r
            max_s = s

            # 当前观察点至目标点与横向格网线的交点数
            x_intersection = math.floor(-x_distance[i] / real_lon)

            # 获取当前目标点单个网格 LOS对应的dx[i] 的dy[j]
            dp_lat = -y_distance[j] / -x_distance[i] * real_lon

            # 格网点对角线长度
            dp_l = math.sqrt(dp_lat ** 2 + real_lon ** 2)

            # 循环遍历x交点，看纵向网格线，从目标点到观察点进行计算的

            k = 1
            while k <= x_intersection:
                # 第k个交点时对应的y长度
                dy = dp_lat * k

                # 交点至观察点线的长度
                length = r - dp_l * k

                lon_index = i + k + min_x

                # 表示目标点在当前纬度格点之间的相对位置
                inner_y = dy % real_lat
                lat_index = j + math.floor(dy / real_lat) + min_y

                # 线性插值
                h = inner_y * (dem.height[lon_index][lat_index + 1] - dem.height[lon_index][lat_index]) / real_lat + \
                    dem.height[lon_index][lat_index] - see_height

                current_s = h / length

                if current_s > max_s:
                    max_s = current_s
                k += 1

            y_intersection = math.floor(-y_distance[j] / real_lat)
            dp_lon = -x_distance[i] / -y_distance[j] * real_lat
            dp_l = math.sqrt(dp_lon ** 2 + real_lat ** 2)
            k = 1
            while k <= y_intersection:
                dx = dp_lon * k
                length = r - dp_l * k

                lat_index = j + k + min_y
                inner_x = dx % real_lon
                lon_index = i + int(dx / real_lon) + min_x

                h = inner_x * (dem.height[lon_index + 1][lat_index] - dem.height[lon_index][lat_index]) / real_lon + \
                    dem.height[lon_index][lat_index] - see_height

                current_s = h / length
                if current_s > max_s:
                    max_s = current_s
                k += 1

            is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                            y_grid_observe, min_x + i, min_y + j, current_quadrant)
            if s >= max_s and is_in_range:
                result[i][j] = 1
            else:
                result[i][j] = 0
            j += 1
        i += 1
    # print(f"第三象限：\n{result}")

    current_quadrant = 4
    # 第四象限
    i = x_grid_count
    while i < width:
        j = 0
        while j < y_grid_count:
            # 观察点到目标点的距离
            r = math.sqrt(x_distance[i] ** 2 + y_distance[j] ** 2)
            # 高度比距离，观察点至目标点仰角
            s = (dem.height[min_x + i][min_y + j] - see_height) / r
            max_s = s
            # 当前观察点至目标点与横向格网线的交点数
            x_intersection = math.floor(x_distance[i] / real_lon)

            # 获取当前目标点单个网格 LOS对应的dx[i] 的dy[j]
            dp_lat = -y_distance[j] / x_distance[i] * real_lon

            # 格网点对角线长度
            dp_l = math.sqrt(dp_lat ** 2 + real_lon ** 2)

            # 循环遍历x交点，看纵向网格线，从目标点到观察点进行计算的
            k = 1
            while k <= x_intersection:
                # 第k个交点时对应的y长度
                dy = dp_lat * k

                # 交点至观察点线的长度
                length = r - dp_l * k

                lon_index = i - k + min_x

                # 表示目标点在当前纬度格点之间的相对位置
                inner_y = dy % real_lat
                lat_index = j + math.floor(dy / real_lat) + min_y

                # 线性插值
                h = inner_y * (dem.height[lon_index][lat_index + 1] - dem.height[lon_index][lat_index]) / real_lat + \
                    dem.height[lon_index][lat_index] - see_height

                current_s = h / length

                if current_s > max_s:
                    max_s = current_s
                k += 1

            y_intersection = math.floor(-y_distance[j] / real_lat)
            dp_lon = x_distance[i] / -y_distance[j] * real_lat
            dp_l = math.sqrt(dp_lon ** 2 + real_lat ** 2)
            k = 1

            while k <= y_intersection:
                dx = dp_lon * k
                length = r - dp_l * k

                lat_index = j + k + min_y
                inner_x = dx % real_lon
                lon_index = i - math.floor(dx / real_lon) + min_x

                h = inner_x * (dem.height[lon_index - 1][lat_index] - dem.height[lon_index][lat_index]) / real_lon + \
                    dem.height[lon_index][lat_index] - see_height

                current_s = h / length
                if current_s > max_s:
                    max_s = current_s
                k += 1

            is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                            y_grid_observe, min_x + i, min_y + j, current_quadrant)
            if s >= max_s and is_in_range:
                result[i][j] = 1
            else:
                result[i][j] = 0
            j += 1
        i += 1
    # print("最终结果：\n", result)

    return result


def analysis_by_xdraw(x_center, y_center, h_center, radius, h_stand, horizontal_start_angle=0.0,
                      horizontal_end_angle=360.0, current_dem=dem):
    global dem
    dem = current_dem
    # 初始化
    x_grid_observe, y_grid_observe, dem_minx, dem_miny, x_grid_count, y_grid_count, see_height, x_grid_center, y_grid_center, min_x, max_x, min_y, max_y, is_in_grid_point, result = get_initial_param_spderl(
        x_center, y_center, h_center, radius, h_stand)
    x_distance, y_distance = get_xy_distance(x_center, y_center, dem_minx, dem_miny, x_grid_count, y_grid_count)

    compare_height = np.zeros((2 * x_grid_count, 2 * y_grid_count))
    # 将观察点四周视为可见
    result[x_grid_count][y_grid_count] = 1
    compare_height[x_grid_count][y_grid_count] = dem.height[x_grid_count + min_x][y_grid_count + min_y] - see_height
    result[x_grid_count - 1][y_grid_count] = 1
    compare_height[x_grid_count - 1][y_grid_count] = dem.height[x_grid_count - 1 + min_x][
                                                         y_grid_count + min_y] - see_height
    result[x_grid_count][y_grid_count - 1] = 1
    compare_height[x_grid_count][y_grid_count - 1] = dem.height[x_grid_count + min_x][
                                                         y_grid_count + min_y - 1] - see_height
    result[x_grid_count - 1][y_grid_count - 1] = 1
    compare_height[x_grid_count - 1][y_grid_count - 1] = dem.height[x_grid_count + min_x - 1][
                                                             y_grid_count + min_y - 1] - see_height

    width = 2 * x_grid_count
    height = 2 * y_grid_count
    real_lon = dem.real_distance
    real_lat = dem.real_distance
    # 当前所在象限
    current_quadrant = 1

    circle_max = x_grid_count
    circle_count = 1
    # 分为8个区
    while circle_count < circle_max:
        # 1:x>0,y>0,    y=rlat,   y++,    x>=y
        current_quadrant = 1
        x_p = x_grid_count + circle_count
        y_p = y_grid_count
        lon_index = x_p + min_x
        while x_p > y_p:
            lat_index = y_p + min_y
            inner_y = y_distance[y_p] / x_distance[x_p] * real_lon
            h_target = dem.height[lon_index][lat_index] - see_height
            # 线性插值 求出上一层对应视线与格网点交点的高度
            h_compare = compare_height[x_p - 1][y_p] - inner_y * (
                    compare_height[x_p - 1][y_p] - compare_height[x_p - 1][y_p - 1]) / real_lat

            # 当前点的距离
            l_target = math.sqrt(y_distance[y_p] ** 2 + x_distance[x_p] ** 2)

            # 上一层点的距离
            l_compare = x_distance[x_p - 1] * l_target / x_distance[x_p]

            # 比较仰角 若上一层仰角大于当前目标点则不可视 同时更新最小可见高程
            if h_compare / l_compare > h_target / l_target:
                compare_height[x_p][y_p] = h_compare / l_compare * l_target
                result[x_p][y_p] = 0
            else:
                compare_height[x_p][y_p] = h_target
                is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                                y_grid_observe, lon_index, lat_index, current_quadrant)
                if is_in_range:
                    result[x_p][y_p] = 1
            y_p += 1

        # 2:x>0,y>0,    x=rlon,   x++,    x<y
        current_quadrant = 1
        x_p = x_grid_count
        y_p = y_grid_count + circle_count
        lat_index = y_p + min_y
        while x_p <= y_p:
            lon_index = x_p + min_x
            inner_x = x_distance[x_p] / y_distance[y_p] * real_lat
            h_target = dem.height[lon_index][lat_index] - see_height
            # 线性插值 求出上一层对应视线与格网点交点的高度
            h_compare = compare_height[x_p][y_p - 1] - inner_x * (
                    compare_height[x_p][y_p - 1] - compare_height[x_p - 1][y_p - 1]) / real_lon

            # 当前点的距离
            l_target = math.sqrt(y_distance[y_p] ** 2 + x_distance[x_p] ** 2)

            # 上一层点的距离
            l_compare = y_distance[y_p - 1] * l_target / y_distance[y_p]

            # 比较仰角 若上一层仰角大于当前目标点则不可视 同时更新最小可见高程
            if h_compare / l_compare > h_target / l_target:
                compare_height[x_p][y_p] = h_compare / l_compare * l_target
                result[x_p][y_p] = 0
            else:
                compare_height[x_p][y_p] = h_target
                is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                                y_grid_observe, lon_index, lat_index, current_quadrant)
                if is_in_range:
                    result[x_p][y_p] = 1
            x_p += 1

        # 3:x<0,y>0,    x=rlon-1, x--,    -x<=y
        current_quadrant = 2
        x_p = x_grid_count - 1
        y_p = y_grid_count + circle_count
        lat_index = y_p + min_y
        while x_p + y_p > 2 * y_grid_count - 1:
            lon_index = x_p + min_x
            inner_x = -x_distance[x_p] / y_distance[y_p] * real_lat
            h_target = dem.height[lon_index][lat_index] - see_height
            # 线性插值 求出上一层对应视线与格网点交点的高度
            h_compare = compare_height[x_p][y_p - 1] - inner_x * (
                    compare_height[x_p][y_p - 1] - compare_height[x_p + 1][y_p - 1]) / real_lon

            # 当前点的距离
            l_target = math.sqrt(y_distance[y_p] ** 2 + x_distance[x_p] ** 2)

            # 上一层点的距离
            l_compare = y_distance[y_p - 1] * l_target / y_distance[y_p]

            # 比较仰角 若上一层仰角大于当前目标点则不可视 同时更新最小可见高程
            if h_compare / l_compare > h_target / l_target:
                compare_height[x_p][y_p] = h_compare / l_compare * l_target
                result[x_p][y_p] = 0
            else:
                compare_height[x_p][y_p] = h_target
                is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                                y_grid_observe, lon_index, lat_index, current_quadrant)
                if is_in_range:
                    result[x_p][y_p] = 1
            x_p -= 1

        # 4:x<0,y>0,    y=rlat,   y++,    -x>y
        current_quadrant = 2
        x_p = (x_grid_count - 1) - circle_count
        y_p = y_grid_count
        lon_index = x_p + min_x
        while x_p + y_p <= 2 * y_grid_count - 1:
            lat_index = y_p + min_y
            inner_y = y_distance[y_p] / -x_distance[x_p] * real_lon
            h_target = dem.height[lon_index][lat_index] - see_height
            # 线性插值 求出上一层对应视线与格网点交点的高度
            h_compare = compare_height[x_p + 1][y_p] - inner_y * (
                    compare_height[x_p + 1][y_p] - compare_height[x_p + 1][y_p - 1]) / real_lat

            # 当前点的距离
            l_target = math.sqrt(y_distance[y_p] ** 2 + x_distance[x_p] ** 2)

            # 上一层点的距离
            l_compare = x_distance[x_p + 1] * l_target / x_distance[x_p]

            # 比较仰角 若上一层仰角大于当前目标点则不可视 同时更新最小可见高程
            if h_compare / l_compare > h_target / l_target:
                compare_height[x_p][y_p] = h_compare / l_compare * l_target
                result[x_p][y_p] = 0
            else:
                compare_height[x_p][y_p] = h_target
                is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                                y_grid_observe, lon_index, lat_index, current_quadrant)
                if is_in_range:
                    result[x_p][y_p] = 1
            y_p += 1

        # 5:x<0,y<0,    y=rlat-1, y--,    -x>=-y
        current_quadrant = 3
        x_p = (x_grid_count - 1) - circle_count
        y_p = y_grid_count - 1
        lon_index = x_p + min_x
        while x_p < y_p:
            lat_index = y_p + min_y
            inner_y = -y_distance[y_p] / -x_distance[x_p] * real_lon
            h_target = dem.height[lon_index][lat_index] - see_height
            # 线性插值 求出上一层对应视线与格网点交点的高度
            h_compare = compare_height[x_p + 1][y_p] - inner_y * (
                    compare_height[x_p + 1][y_p] - compare_height[x_p + 1][y_p + 1]) / real_lat

            # 当前点的距离
            l_target = math.sqrt(y_distance[y_p] ** 2 + x_distance[x_p] ** 2)

            # 上一层点的距离
            l_compare = x_distance[x_p + 1] * l_target / x_distance[x_p]

            # 比较仰角 若上一层仰角大于当前目标点则不可视 同时更新最小可见高程
            if h_compare / l_compare > h_target / l_target:
                compare_height[x_p][y_p] = h_compare / l_compare * l_target
                result[x_p][y_p] = 0
            else:
                compare_height[x_p][y_p] = h_target
                is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                                y_grid_observe, lon_index, lat_index, current_quadrant)
                if is_in_range:
                    result[x_p][y_p] = 1
            y_p -= 1

        # 6:x<0,y<0,    y=rlat-1, y--,    -x>=-y
        current_quadrant = 3
        x_p = x_grid_count - 1
        y_p = (y_grid_count - 1) - circle_count
        lat_index = y_p + min_y
        while x_p >= y_p:
            lon_index = x_p + min_x
            inner_x = -x_distance[x_p] / -y_distance[y_p] * real_lat
            h_target = dem.height[lon_index][lat_index] - see_height
            # 线性插值 求出上一层对应视线与格网点交点的高度
            h_compare = compare_height[x_p][y_p + 1] - inner_x * (
                    compare_height[x_p][y_p + 1] - compare_height[x_p + 1][y_p + 1]) / real_lon

            # 当前点的距离
            l_target = math.sqrt(y_distance[y_p] ** 2 + x_distance[x_p] ** 2)

            # 上一层点的距离
            l_compare = y_distance[y_p + 1] * l_target / y_distance[y_p]

            # 比较仰角 若上一层仰角大于当前目标点则不可视 同时更新最小可见高程
            if h_compare / l_compare > h_target / l_target:
                compare_height[x_p][y_p] = h_compare / l_compare * l_target
                result[x_p][y_p] = 0
            else:
                compare_height[x_p][y_p] = h_target
                is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                                y_grid_observe, lon_index, lat_index, current_quadrant)
                if is_in_range:
                    result[x_p][y_p] = 1
            x_p -= 1

        # 7:x>0,y<0,    x=rlon,   x++,    x<=-y
        current_quadrant = 4
        x_p = x_grid_count
        y_p = (y_grid_count - 1) - circle_count
        lat_index = y_p + min_y
        while x_p + y_p < 2 * y_grid_count - 1:
            lon_index = x_p + min_x
            inner_x = x_distance[x_p] / -y_distance[y_p] * real_lat
            h_target = dem.height[lon_index][lat_index] - see_height
            # 线性插值 求出上一层对应视线与格网点交点的高度
            h_compare = compare_height[x_p][y_p + 1] - inner_x * (
                    compare_height[x_p][y_p + 1] - compare_height[x_p - 1][y_p + 1]) / real_lon

            # 当前点的距离
            l_target = math.sqrt(y_distance[y_p] ** 2 + x_distance[x_p] ** 2)

            # 上一层点的距离
            l_compare = y_distance[y_p + 1] * l_target / y_distance[y_p]

            # 比较仰角 若上一层仰角大于当前目标点则不可视 同时更新最小可见高程
            if h_compare / l_compare > h_target / l_target:
                compare_height[x_p][y_p] = h_compare / l_compare * l_target
                result[x_p][y_p] = 0
            else:
                compare_height[x_p][y_p] = h_target
                is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                                y_grid_observe, lon_index, lat_index, current_quadrant)
                if is_in_range:
                    result[x_p][y_p] = 1
            x_p += 1

        # 8:x>0,y<0,    y=rlat-1, y--,    x>-y
        current_quadrant = 4
        x_p = x_grid_count + circle_count
        y_p = y_grid_count - 1
        lon_index = x_p + min_x
        while x_p + y_p >= 2 * y_grid_count - 1:
            lat_index = y_p + min_y
            inner_y = -y_distance[y_p] / x_distance[x_p] * real_lon
            h_target = dem.height[lon_index][lat_index] - see_height
            # 线性插值 求出上一层对应视线与格网点交点的高度
            h_compare = compare_height[x_p - 1][y_p] - inner_y * (
                    compare_height[x_p - 1][y_p] - compare_height[x_p - 1][y_p + 1]) / real_lat

            # 当前点的距离
            l_target = math.sqrt(y_distance[y_p] ** 2 + x_distance[x_p] ** 2)

            # 上一层点的距离
            l_compare = x_distance[x_p - 1] * l_target / x_distance[x_p]

            # 比较仰角 若上一层仰角大于当前目标点则不可视 同时更新最小可见高程
            if h_compare / l_compare > h_target / l_target:
                compare_height[x_p][y_p] = h_compare / l_compare * l_target
                result[x_p][y_p] = 0
            else:
                compare_height[x_p][y_p] = h_target
                is_in_range = point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_grid_observe,
                                                y_grid_observe, lon_index, lat_index, current_quadrant)
                if is_in_range:
                    result[x_p][y_p] = 1
            y_p -= 1
        circle_count += 1
    # true_result = judge_with_true.get_true_result()
    # print("最终结果是否相同：\n", judge_with_true.are_arrays_equal(result, true_result))
    return result


"""
x_center-观察点横坐标  经度
y_center-观察点纵坐标  纬度
x_grid_center-观察点格网横坐标
y_grid_center-观察点格网纵坐标
h_center-测点地面海拔高程
to_x-边缘点横坐标
to_y-边缘点纵坐标
h_stand-观察者距离地面相对高度
max_x-观察点所处网格位置+半径所需网格点数
"""


def get_initial_param(x_center, y_center, h_center, to_x, to_y, h_stand):
    per_x = dem.dx
    per_y = dem.dy

    # 计算半径
    r = math.sqrt((x_center - to_x) ** 2 + (y_center - to_y) ** 2)
    r_distance = r * dem.rdx
    r_distance2 = r_distance ** 2

    # 某个方向需要的格网点数
    x_grid_count = int(r / dem.dx)
    y_grid_count = int(r / dem.dy)
    see_height = h_center + h_stand

    # 求站立点的格网位置
    x_grid_observe, y_grid_observe, x_grid_center, y_grid_center, x_grid_inner, y_grid_inner = dem.get_point_location(
        x_center, y_center)

    min_x: int = x_grid_center - (x_grid_count - 1)
    max_x: int = x_grid_center + x_grid_count
    min_y: int = y_grid_center - (y_grid_count - 1)
    max_y: int = y_grid_center + y_grid_count
    if min_x < 0 or min_y < 0 or max_x >= dem.x_size or max_y >= dem.y_size:
        print("数据范围不在求解区域")

    # 第一个格网点的经纬度
    dem_minx = x_center - x_grid_inner - (x_grid_count - 1) * dem.dx
    dem_miny = y_center - y_grid_inner - (y_grid_count - 1) * dem.dy

    # result数组是以中心网格点为中心建立的二维数组 不一定是中心观察点
    global result

    result = np.zeros((2 * x_grid_count, 2 * y_grid_count))

    return dem_minx, dem_miny, per_x, per_y, x_grid_count, y_grid_count, see_height, x_grid_center, y_grid_center, min_x, max_x, min_y, max_y


# spderl算法变量初始化
def get_initial_param_spderl(x_center, y_center, h_center, radius, h_stand):
    # 实际半径radius 转为 与格网间距同意尺度半径
    r_grid = radius / dem.rdx

    # 某个方向需要的格网点数
    x_grid_count = int(r_grid / dem.dx)
    y_grid_count = int(r_grid / dem.dy)
    see_height = h_center + h_stand

    # 求站立点的格网位置
    x_grid_observe, y_grid_observe, x_grid_center, y_grid_center, x_grid_inner, y_grid_inner, is_in_grid_point = dem.get_point_location(
        x_center, y_center)

    # 第一个格网点的经纬度
    dem_minx = x_center - x_grid_inner - (x_grid_count - 1) * dem.dx
    dem_miny = y_center - y_grid_inner - (y_grid_count - 1) * dem.dy

    min_x: int = x_grid_center - (x_grid_count - 1)
    max_x: int = x_grid_center + x_grid_count
    min_y: int = y_grid_center - (y_grid_count - 1)
    max_y: int = y_grid_center + y_grid_count
    if min_x < 0 or min_y < 0 or max_x >= dem.x_size or max_y >= dem.y_size:
        print("数据范围不在求解区域")

    # result数组是以中心网格点为中心建立的二维数组 不一定是中心观察点
    # global result

    result = np.zeros((2 * x_grid_count, 2 * y_grid_count))

    return x_grid_observe, y_grid_observe, dem_minx, dem_miny, x_grid_count, y_grid_count, see_height, x_grid_center, y_grid_center, min_x, max_x, min_y, max_y, is_in_grid_point, result


def get_start_loc(x_center, y_center, radius):
    # 实际半径radius 转为 与格网间距同意尺度半径
    r_grid = radius / dem.rdx

    # 某个方向需要的格网点数
    x_grid_count = int(r_grid / dem.dx)
    y_grid_count = int(r_grid / dem.dy)

    # 求站立点的格网位置
    x_grid_observe, y_grid_observe, x_grid_center, y_grid_center, x_grid_inner, y_grid_inner, _ = dem.get_point_location(
        x_center, y_center)

    # 第一个格网点的经纬度-右上角
    start_lon = x_center - x_grid_inner - (x_grid_count - 1) * dem.dx
    start_lat = y_center + y_grid_inner + (y_grid_count - 1) * dem.dy
    return start_lon, start_lat


# 获取中心点距离各格网线的距离数组(在DEM定义的投影中)
# 这里x_center为实际观察者的经纬度值 因此对应的p和d均为实际值 而不是格网点值
def get_xy_distance(x_center, y_center, dem_min_x, dem_min_y, x_grid_count, y_grid_count):
    x_distance = np.zeros(2 * x_grid_count)
    y_distance = np.zeros(2 * y_grid_count)
    for i in range(2 * x_grid_count):
        x_distance[i] = (dem_min_x - x_center + i * dem.dx) * dem.rdx
    for j in range(2 * y_grid_count):
        y_distance[j] = (dem_min_y - y_center + j * dem.dy) * dem.rdy
    return x_distance, y_distance


# 根据角度值及观察点对应格网中位置生成对应的边界直线方程
def calculate_line_equation(angle, x, y):
    quadrant = judge_angle_at_quadrant(angle)
    # 根据输入角度转为对应区间与x夹角值 ×会出现错误
    # 输入角度转为对应与 x 正方向夹角值
    if quadrant == 1:
        angle = 90 - angle
    elif quadrant == 2:
        angle = 450 - angle
    elif quadrant == 3:
        angle = 270 - angle
    else:
        angle = 270 - angle
    # 将角度转换为弧度
    if angle == 90 or angle == 270:
        return -1, -1
    angle_rad = math.radians(angle)

    # 计算斜率
    slope = math.tan(angle_rad)
    # 设置一个阈值，如果计算结果接近于0，则视为0
    threshold = 1e-10  # 例如, 1e-10
    if abs(slope) < threshold:
        slope = 0

    # 解方程 y = mx + b，计算截距 b
    b = y - slope * x

    # 返回斜率和截距
    return slope, b


# 根据当前格网索引(非结果数组)和角度起止范围判断是否处于可视距内,需要结合象限进行考虑
def point_is_in_range(horizontal_start_angle, horizontal_end_angle, x_observe, y_observe, x, y, current_quadrant):
    angle = 0
    # 分象限考虑 先计算当前坐标与靠近北方向的夹角

    # 位于1、3象限 需计算与y轴夹角
    if current_quadrant == 1 or current_quadrant == 3:
        # 返回弧度制
        # 考虑同y轴的情况
        if y != y_observe:
            atan_value = math.atan(abs((x - x_observe) / (y - y_observe))) + quadrant[current_quadrant - 1]
            angle = math.degrees(atan_value) % 360
        else:
            angle = 90 if x > x_observe else 270
    # 转为角度值 可能存在负数的情况

    # 位于2、4象限 需计算与x夹角
    elif current_quadrant == 2 or current_quadrant == 4:
        if x != x_observe:
            atan_value = math.atan(abs((y - y_observe) / (x - x_observe))) + quadrant[current_quadrant - 1]
            angle = math.degrees(atan_value) % 360
        else:
            angle = 0 if y > x_observe else 180
    if horizontal_start_angle <= angle <= horizontal_end_angle:
        return True
    else:
        return False


# 根据角度返回对应笛卡尔坐标系的象限
def judge_angle_at_quadrant(angle):
    if 0 <= angle < 90:
        quadrant = 1
    elif 90 <= angle < 180:
        quadrant = 4
    elif 180 <= angle < 270:
        quadrant = 3
    else:
        quadrant = 2
    return quadrant


def judge_is_out_bound(x_center, y_center, radius, current_dem):
    is_out_bound = False

    # 实际半径radius 转为 与格网间距同意尺度半径
    r_grid = radius / current_dem.rdx

    # 某个方向需要的格网点数
    x_grid_count = int(r_grid / current_dem.dx)
    y_grid_count = int(r_grid / current_dem.dy)

    # 求站立点的格网位置
    x_grid_observe, y_grid_observe, x_grid_center, y_grid_center, x_grid_inner, y_grid_inner, is_in_grid_point = current_dem.get_point_location(
        x_center, y_center)

    min_x: int = x_grid_center - (x_grid_count - 1)
    max_x: int = x_grid_center + x_grid_count
    min_y: int = y_grid_center - (y_grid_count - 1)
    max_y: int = y_grid_center + y_grid_count
    if min_x < 0 or min_y < 0 or max_x >= current_dem.x_size or max_y >= current_dem.y_size:
        is_out_bound = True
    return is_out_bound


# 当前模块是否作为主程序直接运行
if __name__ == '__main__':
    # main函数中声明的变量默认为global variable
    # 而其他def函数中声明的变量则默认为local variable

    # dem文件路径
    file_path = r"C:\Users\GuoQing\Desktop\Experiment\S_PDERL\data\Copernicus_DSM_COG_10_N28_00_E097_00_DEM.tif"
    # 分辨率
    real_distance = 30
    # 实例化dem类
    dem = dem_data.Dem(file_path, real_distance)

    # 保存对应象限需要添加的角度
    quadrant = [0, math.pi * 3 / 2, math.pi, math.pi / 2]

    # 结果矩阵 初始化为0
    # 其中 result[x,y]表示的是第x行y列 对应在XY坐标系中是x->纵轴 y->横纵
    result = []

    x_center = 97.4157
    y_center = 28.711
    h_center = 3978.34673937845
    to_x = 97.421364
    to_y = 28.706327
    h_stand = 5
    radius = 793.030586294823

    horizontal_start_angle = 0
    horizontal_end_angle = 360
    # print("PDERL:\n")
    # # analysis_by_pderl(x_center, y_center, h_center, to_x, to_y, h_stand)

    print("XPDERL:\n")
    analysis_by_xpderl_with_angle(97.86274697268131, 28.863024750459086, 4000, 0.1 * dem.rdx, 1, 359, 150)
    #
    # print("R3:\n")
    # analysis_by_r3(x_center, y_center, h_center, radius, h_stand, horizontal_start_angle, horizontal_end_angle)
    #
    # print("SPDERL:\n")
    # analysis_by_spderl(x_center, y_center, h_center, radius, h_stand, horizontal_start_angle, horizontal_end_angle)
    # print("Xdraw:\n")
    # xdraw_res = analysis_by_xdraw(x_center, y_center, h_center, radius, h_stand, horizontal_start_angle,
    #                               horizontal_end_angle, dem)
    # print(xdraw_res)

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
    print("右面：\n", result)

    true_result = judge_with_true.get_true_result()
    print("右面是否相同：\n", judge_with_true.are_arrays_equal(result, true_result))
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
    print("上面：\n", result)

    # 比较与原算法结果是否有异
    # true_result = judge_with_true.get_true_result()
    # print("是否相同：\n", judge_with_true.are_arrays_equal(result, true_result))

    true_result = judge_with_true.get_true_result()
    print("上面是否相同：\n", judge_with_true.are_arrays_equal(result, true_result))

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
    print("左面：\n", result)

    # 比较与原算法结果是否有异
    true_result = judge_with_true.get_true_result()
    print("左面是否相同：\n", judge_with_true.are_arrays_equal(result, true_result))
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

    print("最终结果:\n", result)

    # 比较与原算法结果是否有异
    true_result = judge_with_true.get_true_result()
    print("最终是否相同：\n", judge_with_true.are_arrays_equal(result, true_result))
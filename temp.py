def analysis_by_r3(x_center, y_center, h_center, radius, h_stand, horizontal_start_angle=0.0,
                   horizontal_end_angle=360.0, current_dem=dem):
    global dem
    dem = current_dem
    # 初始化
    x_grid_observe, y_grid_observe, dem_minx, dem_miny, x_grid_count, y_grid_count, see_height, x_grid_center, y_grid_center, min_x, max_x, min_y, max_y = get_initial_param_spderl(
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

    # 将不再分析区域内的值转为0

    true_result = judge_with_true.get_true_result()
    # print("最终结果是否相同：\n", judge_with_true.are_arrays_equal(result, true_result))
    return result
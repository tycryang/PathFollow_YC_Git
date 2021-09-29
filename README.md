自我介绍：车辆工程专业毕业，主要方向是车辆动力学仿真与控制，操纵稳定性范畴，对操纵动力学、悬架K&C特性、轮胎F&M特性比较了解吧，毕业后在主机厂做性能集成，感觉部门实力太弱，想换方向，思前想后想去觉得不如往风口上靠——做自动驾驶吧，同时为了发挥自己的优势，决定转行方向是自动驾驶控制算法工程师，所以便有了这么一个项目。

项目简介：（1）使用PID、LQR、MPC算法实现对各种轨迹的跟踪；同时车辆闭环仿真有一个驾驶员预瞄模型，也可以用于做轨迹跟踪，也一并实现了。分成了两个文件夹：PID_and_DriverModel，实现了PID和DriverModel利用侧向偏差实现轨迹跟踪，LQRnMPC_LookAheadorNot，里面有LQR、MPC、采用增量输入的MPC和考虑软约束的MPC四个控制器。
为什么会有后面两个控制器？LQR和MPC是最先实现的，实现后发现轨迹可以跟的大差不差，但是方向盘转角波动很大，肯定影响工程化实现。后根据北理龚建伟老师的《无人驾驶车辆模型预测仿真控制》将控制量从转角变为转角增量，这样就可以增加对转角速度的约束，取得了一些进步。同时《无人驾驶车辆模型预测仿真控制》还给了一个加软约束的建议，虽然此书没有写的很清楚，但是经过一番探索还是做出来了，这个软约束其实对控制而言并不是一个有利的因素，这是一个为了提高乘客舒适性的优化方向，比方说车辆侧向加速度最高可达0.9g，但是日常没谁用，日常通常不会超过0.4g，但是你不可以直接把车辆限制在0.4g以内，万一需要紧急壁障呢？我觉得这个软约束就是用来解决这个矛盾的，超出0.4g的侧向加速度写成代价函数。总体实现效果跟这些预期差不多。对转角增量增加约束后方向盘转角波动会变小很多，增加软约束后方向盘转角的突变变少了/轨迹变更平滑了（尤其紧急避障这种）

但是转角依然不够平滑，所以我又想到了驾驶员预瞄模型，这个涉及偏差计算环节，网上很多apollo普及文章告诉我们，偏差是基于最近点计算的（更细一点就是比最近点更近的投影点），或者再说的直白一点，反馈控制是有了偏差才去控制，要控制自然要变化控制量；哪怕考虑了路径曲率，但是也是一个点上的曲率（或许也考虑用前方一段路径做一个平均曲率也可以优化转角平滑度，所以我也尝试对规划轨迹的曲率进行smooth，有一些效果吧）。而来自车汽车闭环仿真的驾驶员预瞄跟随模型，原理是模拟驾驶员开车视野看的车前方一段距离，我是希望在未来一个较短时间内达到前方想要的轨迹上，而不是现在，因为乘用车也没有单独控制侧向平动的输入。而且郭孔辉院士的《操纵动力学原理》中也举了一个很典型的例子帮助我们理解这个驾驶员模型中预瞄的意义，新手驾驶员开车紧张，生怕产生一点点偏差，盯着当前偏差对方向进行修正，结果往往在车道内画龙；而老司机通常视野更远，操作方向盘的修正操作也更少。基于此在LQRnMPC_LookAheadorNot中设计了两种偏差计算模型，一种当然是基于当前偏差，另一种是引入预瞄时间，基于未来0.1~0.2s后的偏差；此时方向盘转角的波动又变小很多了，总体我觉得应该还可以。

以上，当然里面还有些东西可能理解的还不够到位，比方说其实对规划轨迹曲率smooth后，而且定曲率了为什么还有较大的抖动，抖动来自哪里其实我还是没有完全想通？在比方说，方向角及偏差范围限定在-pi,pi之间，如果其中一个从pi突变到-pi而另一个没有跟上，此时产生的偏差要怎么处理，有什么好的办法？如果有幸你有想法，欢迎交流。邮箱：yangchao9103@163.com

所使用的工具：matlab 2019b，编程语言：m语言，Simulink建模自建模块尽量使用的sfunction。

后续优化方向：
（1）将现有的车辆模型从线性2DOF模型加上非线性的轮胎模型，抑或换成更多自由度14DOF模型，甚至CarSim模型，尤其是CarSim模型还可以表达一些车辆侧偏动力学中的延迟，可以进一步发挥预瞄的作用。
（2）在车辆动力学的输出中加入一些扰动，模拟实际传感器采集信号有噪声的情况，然后利用低通滤波/卡尔曼滤波去进行滤波，控制效果又会面临什么新的课题呢？

最后，愿自己可以成功换方向！

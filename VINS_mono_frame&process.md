# VINS_mono



## 整体框架

![image-20210825073742749](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210825073742749.png)

1、前端(数据预处理)

图像：提取图像 Harris 角点,利用金字塔光流跟踪相邻帧,通过 RANSAC 去除异常点,最后将跟踪到的特征点 push 到图像队列中,并通知后端进行处理。

IMU：将 IMU 数据进行积分,得到当前时刻的位置、速度和旋转(PVQ),同时计算在后端优化中将用到的相邻帧的预积分增量,及预积分误差的 Jacobian 矩阵和协方差项。

相邻帧的光流跟踪：建立当前帧与滑窗中其他帧的共视关系,即重投影误差约束

Shi-Tomas特征点提取：保证一定数量的特征点

IMU预积分：计算当前帧的位姿作为初始值,计算IMU约束的误差项、协方差和Jacobian

2、后端(滑窗优化)

将视觉约束、IMU 约束和闭环约束放在一个大的目标函数中进行非线性优化,求解滑窗内所有帧的 PVQ、bias 等。

联合优化视觉重投影误差、IMU运动约束、先验约束

3、初始化

首先,利用 SFM 进行纯视觉估计滑窗内所有帧的位姿及 3D 点逆深度,最后与 IMU 预积分进行对齐求解初始化参数。

基于视觉信息来计算运动轨迹

基于IMU信息来计算运动轨迹

视觉和IMU松耦合来优化尺度、速度、重力

4、闭环

利用 DBoW 进行闭环检测,当检测成功后进行重定位,最后对整个相机轨迹进行闭环优化。

通过词袋找闭环帧,通过暴力匹配进行特征点匹配,通过滑窗优化闭环帧的相对位姿

利用闭环帧的相对位姿进行轨迹的四自由度优化



## ROS节点

![image-20210825074325263](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210825074325263.png)

```
$ rosrun rqt_graph rqt_graph
```



## ROS基础知识

master——节点管理器,负责管理所有节点

node——节点,相互独立的可执行文件

topic——主题,节点之间相互通信的主题



rosrun、roslaunch、package、node

rosrun: 一次运行一个node
roslaunch: 一次运行多个node
package: 一个功能包,包含多个node
node: 一个节点就是package中的一个可执行文件

运行launch文件,以.launch为后缀名,放在Package的launch目录内。

```
$ roslaunch [PackageName] [LaunchFileName]
```

如何运行VINS-Mono:

```
$ roslaunch vins_estimator euroc.launch
$ roslaunch vins_estimator vins_rviz.launch
$ rosbag play YOUR_PATH_TO_DATASET/MH_01_easy.bag
```

### euroc.launch: 为了一次运行多个node

![image-20210825074701717](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210825074701717.png)

launch文件解释：

pkg=package名字(在package.xml定义)

type=可执行文件名字(在CMakeLists.txt定义)

name=node名字，相当于ros::init()



可以看到包含的节点有：

feature_tracker、vins_estimator、pose_graph

node1: 前端

node2: 后端

node3: 闭环



### package.xml 功能包:

![image-20210825075000459](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210825075000459.png)

<name> - 名字
<version> - 版本
<description> - 内容描述
<maintainer> - 作者
<license> - 软件发行版通行证

<buildtool_depend> - 指定编译工具
<build_depend> - 指定头文件、链接库、其他源文件
<run_depend> - 指定运行所需的其他功能包
<test_depend> - 指定单元测试所需的其他功能包



### CMakeLists.txt:学习一个工程,从CMakeLists开始

![image-20210825075142489](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210825075142489.png)

project名字需与package相同



![image-20210825075332320](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210825075332320.png)



### 配置文件:euroc_config.yaml

最大的点数: max_cnt

最大的求解时间: max_solver_time

最大的迭代此书: max_num_iterations

IMU噪声:

```
acc_n: 0.2          # accelerometer measurement noise standard deviation. #0.2
gyr_n: 0.02         # gyroscope measurement noise standard deviation.     #0.05
acc_w: 0.0002         # accelerometer bias random work noise standard deviation.  #0.02
gyr_w: 2.0e-5       # gyroscope bias random work noise standard deviation.     #4.0e-5
g_norm: 9.81007     # gravity magnitude
```

IMU和Cam的时间戳同步:

```
estimate_td: 0                      # online estimate time offset between camera and imu
td: 0.0                             # initial value of time offset. unit: s. readed image clock + td = real image clock (IMU clock)
```



## 后台estimator_node.cpp



![image-20210825080451401](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210825080451401.png)

### 后台IMU、图像、后端三个线程，线程互斥锁

<img src="%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210825081001576.png" alt="image-20210825081001576" style="zoom:80%;" />

<img src="%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210825081010973.png" alt="image-20210825081010973" style="zoom:80%;" />

mutex 多线程互斥锁 m_buf

condition_variable 条件变量 con

• IMU、图像、后端三个线程共用一个m_buf

• 后端process()线程调用wait(),自动调用m_buf.lock()来加锁;

• 若getMeasurements未得到图像和两帧之间的IMU,则返回false,线程被阻塞,此时,wait()会自动调用m_buf.unlock()来释放锁,使得imu和feature的回调线程继续push;

• imu和feature的回调线程,每拿到数据,都会调用con.notify_one唤醒process()线程,使其继续尝试
getMeasurements。



## 前端Node: feature_tracker

前端main():启动前端回调 img_callback()

### readImage()

![image-20210825081609270](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210825081609270.png)

![image-20210825081754404](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210825081754404.png)

cur_img 上一帧图像
forw_img 当前帧图像
cur_pts 上一帧的点坐标
forw_pts 当前帧的点坐标
ids 每个点的id号
track_cnt 每个点被跟踪的次数
cur_un_pts 最新帧的归一化相机系坐标

### 前端发布出去的消息

![image-20210830075433950](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210830075433950.png)

topic: feature

![image-20210830075501514](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210830075501514.png)

```
[sensor_msgs/PointCloud]:
❑ header
时间戳
-seq
-stamp
-frame_id
❑ points 归一化相机系坐标
-x、y、z
❑ channels
(id、u、v、v x 、v y )
-name
-values
```



## 后端Node: vins_estimator



### getMeasurements

根据时间戳,挑选当前帧和上一帧图像之间的imu数据,用于后续的IMU积分

![image-20210830080042417](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210830080042417.png)



挑选比第一张图像老的imu数据

![image-20210830080157930](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210830080157930.png)

上述3种情形：

1 若IMU太老,则等一下IMU；

2 若图像太老,则删除最老的图像；

3 若不新不旧,则携手双双把家还。

则两个IMU信息和一帧图像组合，用来进行预积分和优化：

![image-20210830080310748](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210830080310748.png)



### 视觉和IMU的联合优化

![image-20210830080637690](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210830080637690.png)

IMU预积分作用:计算IMU观测的误差、协方差、Jacobian

协方差的根源是因为imu有噪声,因此随着预测的增加,观测的不确定度会越来越大,即协方差越來越大;
当来了视觉的观测后,系统会根据两者的协方差大小,来权衡更相信谁。

![image-20210830080732994](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210830080732994.png)

### processIMU: IMU预积分

拿到了图像和两帧之间的多个imu信息



```
pre_integrations[frame_count]->push_back(dt, linear_acceleration, angular_velocity);
```

push_back进行预积分,根据两帧之间的IMU信息,计算两帧之间的位置、旋转、速度的变化量,作为IMU约束的误差项,并根据IMU的噪声大小,计算IMU约束的协方差(即不确定度)



```
		int j = frame_count;         
        Vector3d un_acc_0 = Rs[j] * (acc_0 - Bas[j]) - g;
        Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - Bgs[j];
        Rs[j] *= Utility::deltaQ(un_gyr * dt).toRotationMatrix();
        Vector3d un_acc_1 = Rs[j] * (linear_acceleration - Bas[j]) - g;
        Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
        Ps[j] += dt * Vs[j] + 0.5 * dt * dt * un_acc;
        Vs[j] += dt * un_acc;
```

中值积分预测最新帧的位姿,作为视觉的初始位姿,根据初始位姿,重投影，计算重投影误差



### IntegrationBase::push_back() IMU预积分

![image-20210830081925640](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210830081925640.png)

push_back(dt, acc, gyr)

propagate

midPointIntegration 更新IMU约束的误差项、协方差、Jacobian



### midPointIntegration

更新IMU约束的误差项

![image-20210830082051651](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210830082051651.png)

![image-20210830082040898](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210830082040898.png)



更新IMU约束的协方差、Jacobian

![image-20210830082125945](%E5%9B%BE%E8%A1%A8%E5%BA%93/image-20210830082125945.png)
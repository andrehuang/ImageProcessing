关于代码的说明：

1 deblurring问题中模糊算子的边值条件为'circular'；

2 Perona-Malik PDE实现中的c(s)实现了两种：1/(1+s/K) 和 1/(1+(s/K)^2)，可以用fchoice来选择；

3 shock filter问题中，实现的F(s) = sign(s)；L(u)实现了要求的两种，可以用Lchoice来选择。

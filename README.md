# Лабораторные работы по дисциплине "Методы оптимизации" на факультете ФПМИ, НГТУ


### Папки:
### 1. Методы одномерного поиска
> Реализовать методы дихотомии, золотого сечения, исследовать их сходимость и провести сравнение по числу вычислений функции для 
достижения заданной точности **eps** от **10^(-1)** до **10^(-7)**.  
>Построить график зависимости количества вычислений минимизируемой функции от десятичного логарифма задаваемой точности.  
>Реализовать алгоритм поиска интервала, содержащего минимум функции. Реализовать метод Фибоначчи, сравнить его с методами дихотомии и 
золотого сечения.  
>Задачу выполнять для функции: **f(x) = (x - 2)^2, x ϵ [-2, 20]**.

### 2. Методы спуска (0-го, 1-го и 2-го порядка и переменной метрики)
> Реализовать два метода поиска экстремума функции: **метод Розенброка**, **метод Бройдена**. Включить в реализуемый алгоритм собственную процедуру, реализующую одномерный поиск по направлению.  
> С использованием разработанного программного обеспечения исследовать алгоритмы на квадратичной функции
> <p align="center"> <img width="440" src="2-quadratic-function.png"> </p>
> , функции Розенброка 
> <p align="center"> <img width="440" src="2-rosenbrock-function.png"> </p>
> и на заданной в соответствии с вариантом тестовой функции, осуществляя спуск из различных исходных точек (не менее двух). Исследовать сходимость алгоритма, фиксируя точность определения минимума/максимума, количество итераций метода и количество вычислений функции в зависимости от задаваемой точности поиска. Результатом выполнения данного пункта должны быть выводы об объёме вычислений в зависимости от задаваемой точности и начального приближения.

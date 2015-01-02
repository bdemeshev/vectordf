vectordf
========

Пакет `vectordf` предназначен для работы с распределениями, описываемыми вектором степеней свободы (t, гамма).

Установить пакет можно с помощью команд:

```r
install.packages("devtools")
library("devtools")
install_github("bdemeshev/vectordf")
```


Примеры использования.

Матричное нормальное распределение:
```r
library("sophisthse")
X <- rmatnorm(n = 5, Mu = matrix(c(0,0,1,-1), nrow=2))
```



Планы: 
* дописать функции плотности
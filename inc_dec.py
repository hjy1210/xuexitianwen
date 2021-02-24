def increment(x):
    return x + 1

def decrement(x):
    return x - 1

def fib(n):
    if not type(n) is int:
        raise TypeError("should be an integer")
    if n <=0:
        raise ValueError("should be greater than zero")
    if n==1:
        return 1
    elif n==2:
        return 1
    else:
        return fib(n-1)+fib(n-2)


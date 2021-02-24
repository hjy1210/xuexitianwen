import inc_dec    # The code to test
import pytest

def test_increment():
    assert inc_dec.increment(3) == 4

def test_decrement():
    assert inc_dec.decrement(3) == 2

def test_fib():
    assert inc_dec.fib(3)==2
    with pytest.raises(ValueError):
        inc_dec.fib(-1)
    with pytest.raises(TypeError):
        inc_dec.fib(3.1)

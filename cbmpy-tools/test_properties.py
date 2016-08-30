import weakref

class Reagent():
    _value = None
    _weakref_ = False
    
    @property
    def coefficient(self):
        print(self._value, self._weakref_)

    @coefficient.setter
    def coefficient(self, value):
        if type(value) is Parameter:
            self._value = weakref.ref(value)
            self._weakref_ = True
        else:
            self._value = value
            self._weakref_ = False
    
    @coefficient.getter
    def coefficient(self):
        if self._weakref_:
            x = self._value().value
        else:
            x = self._value
        return x

class Parameter(object):
    _value = None

    #def __call__(self):
        #return self._value
    
    @property
    def value(self):
        pass

    @value.setter    
    def value(self, value):
        self._value = value
 
    @value.getter    
    def value(self):
        return self._value

if __name__ == '__main__':
    p1 = Parameter()
    p2 = Parameter()
    p1.value = 3
    p2.value = 4
    
    r1 = Reagent()
    r2 = Reagent()
    r3 = Reagent()
    r1.coefficient = 5
    r2.coefficient = p1
    r3.coefficient = p2
    

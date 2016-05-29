'''
interval.py - Deal with intervals.
author: Xiao-Ou Zhang <kepbod@gmail.com>
version: 1.0
'''

import copy


class Interval(object):
    '''
    Class: Interval

    Maintainer: Xiao-Ou Zhang

    Version: 1.0

    Usage: a = Interval(list)
           (nested list: [[x,x,f1...],[x,x,f2...]...] / [[x,x],[x,x]...] or
            simple list: [x,x,f1...] / [x,x])
    Notes: all the intervals in the list will become mutually exclusive and
           be sorted after instantiation.

    For example: input:  [[1, 10, 'a'], [17, 22, 'b'], [7, 12, 'c'], [20, 25, 'd'], [30, 35, 'e']]
                 output: [[1, 12, 'a', 'c'], [17, 25, 'b', 'd'], [30, 35, 'e']]

    Attributes: interval

    Functions: c = a + b or a += b
               c = b + a
               c = a * b or a *= b
               c = b * a
               c = a - b or a -= b
               c = b - a
               a[n] or a[n:m]
               [x, x] in a or [[x, x], [x, x]] not in a
               a.complement(sta, end)
               a.extractwith(b)
               a.extractwithout(b)
               mapto(interval, index) -> interval
               overlapwith(index, interval) -> index
    '''
    def __init__(self, interval, instance_flag=0):
        self.interval = [[int(i[0]), int(i[1])] + i[2:]
                         for i in Interval.__convert(interval)]
        if not self.interval:
            return
        if not instance_flag:
            self.interval.sort()
            tmp = []
            a = self.interval[0]
            for b in self.interval[1:]:
                if a[1] <= b[0]:
                    tmp.append(a)
                    a = b
                else:
                    a[1] = b[1] if b[1] > a[1] else a[1]
                    a.extend(b[2:])
            tmp.append(a)
            self.interval = tmp

    def __add__(self, interval):
        '''
        Usage: c = a + b or a += b
        extract union intervals, 'a' should be instance.
        '''
        tmp = copy.deepcopy(self.interval)
        if isinstance(interval, Interval):
            tmp.extend(interval.interval)
        else:
            tmp.extend(Interval.__convert(interval))
        return Interval(tmp)

    def __radd__(self, interval):
        '''
        Usage: c = b + a
        extract union intervals, 'a' should be instance.
        '''
        return self.__add__(interval)

    def __mul__(self, interval, real_flag=1):
        '''
        Usage: c = a * b or a *= b
        extract intersection intervals, 'a' should be instance.
        '''
        tmp = []
        tmp1 = self.interval
        if isinstance(interval, Interval):
            tmp2 = interval.interval
        else:
            tmp2 = Interval(interval).interval
        if not tmp1 or not tmp2:
            return Interval([])
        a, b = tmp1[0], tmp2[0]
        i, j = 1, 1
        while True:
            sta = a[0] if a[0] > b[0] else b[0]
            end = a[1] if a[1] < b[1] else b[1]
            if sta < end:
                if real_flag:
                    tmp.append([sta, end] + a[2:] + b[2:])
                else:
                    tmp.append(copy.copy(a))
            if a[1] == end:
                if i == len(tmp1):
                    break
                a = tmp1[i]
                i += 1
            if b[1] == end:
                if j == len(tmp2):
                    break
                b = tmp2[j]
                j += 1
        return Interval(tmp, 1)

    def __rmul__(self, interval):
        '''
        Usage: c = b * a
        extract intersection intervals, 'a' should be instance.
        '''
        return self.__mul__(interval)

    def __sub__(self, interval, real_flag=1):
        '''
        Usage: c = a - b or a -= b
        extract difference intervals, 'a' should be instance.
        '''
        if not self.interval:
            return Interval([])
        if isinstance(interval, Interval):
            tmp = copy.deepcopy(interval)
        else:
            tmp = Interval(interval)
        if not tmp:
            return copy.deepcopy(self)
        if self.interval[0][0] < tmp.interval[0][0]:
            sta = self.interval[0][0]
        else:
            sta = tmp.interval[0][0]
        if self.interval[-1][1] > tmp.interval[-1][1]:
            end = self.interval[-1][1]
        else:
            end = tmp.interval[-1][1]
        tmp.complement(sta, end)
        return self.__mul__(tmp, real_flag)

    def __rsub__(self, interval):
        '''
        Usage: c = b - a
        extract difference intervals, 'a' should be instance.
        '''
        if isinstance(interval, Interval):
            tmp = copy.deepcopy(interval)
        else:
            tmp = Interval(interval)
        if not self.interval:
            return tmp
        if not tmp:
            return Interval([])
        if self.interval[0][0] < tmp.interval[0][0]:
            sta = self.interval[0][0]
        else:
            sta = tmp.interval[0][0]
        if self.interval[-1][1] > tmp.interval[-1][1]:
            end = self.interval[-1][1]
        else:
            end = tmp.interval[-1][1]
        tmp_a = copy.deepcopy(self)
        tmp_a.complement(sta, end)
        return Interval.__mul__(tmp, tmp_a)

    def __getitem__(self, index):
        '''
        Usage: a[n] or a[n:m]
        intercept index and slice on interval objects.
        '''
        return self.interval[index]

    def __repr__(self):
        '''
        print objects.
        '''
        return repr(self.interval)

    def __contains__(self, interval):
        '''
        Usage: [x, x] in a or [[x, x], [x, x]] not in a
        judge whether interval is in a or not, 'a' should be instance.
        '''
        tmp = self.__mul__(interval).interval
        if tmp:
            return True
        else:
            return False

    def complement(self, sta='#', end='#'):
        '''
        Usage: a.complement(sta, end)
        complement of 'a'.
        '''
        tmp = []
        if sta != '#' and sta < self.interval[0][0]:
            tmp.append([sta, self.interval[0][0]])
        a = self.interval[0][1]
        for item in self.interval[1:]:
            b = item[0]
            if a != b:
                tmp.append([a, b])
            a = item[1]
        if end != '#' and end > a:
            tmp.append([a, end])
        self.interval = tmp

    def extractwith(self, interval):
        '''
        Usage: a.extractwith(b)
        extract intervals in 'b'.
        '''
        self.interval = self.__mul__(interval, 0).interval

    def extractwithout(self, interval):
        '''
        Usage: a.extractwithout(b)
        extract intervals not in 'b'.
        '''
        self.interval = self.__sub__(interval, 0).interval

    @staticmethod
    def mapto(interval, index):
        '''
        mapto(interval, index) -> interval
        Map interval onto index.
        '''
        tmp1 = Interval.__init(interval)
        tmp2 = Interval.__init(index)
        return Interval.__map(tmp2, tmp1, flag=1)

    @staticmethod
    def overlapwith(index, interval):
        '''
        overlapwith(index, interval) -> index
        Overlap index with interval.
        '''
        tmp1 = Interval.__init(index)
        tmp2 = Interval.__init(interval)
        return Interval.__map(tmp1, tmp2, flag=0)

    @staticmethod
    def __convert(interval):
        assert type(interval) is list, 'the type you use is {}'.format(type(interval))
        if not interval:
            return interval
        if type(interval[0]) is list:
            return interval
        else:
            return [interval]

    @staticmethod
    def __init(interval):
        mapping = [[int(i[0]), int(i[1])] + i[2:] for i in Interval.__convert(interval)]
        mapping.sort()
        return mapping

    @staticmethod
    def __map(index, interval, flag):
        mapped_fragment = []
        tmp_fragment = []
        if not interval:
            if flag:
                return mapped_fragment
            else:
                return index
        for dex in index:
            dex_info = dex[2:]
            while True:
                try:
                    fragment = interval.pop(0)
                except IndexError:
                    if tmp_fragment:
                        interval.extend(tmp_fragment)
                        tmp_fragment = []
                        continue
                    else:
                        if flag:
                            return mapped_fragment
                        else:
                            return index
                if fragment[0] >= dex[1]:
                    interval.insert(0, fragment)
                    interval[0:0] = tmp_fragment
                    tmp_fragment = []
                    break
                elif dex[0] < fragment[1] and dex[1] > fragment[0]:
                    dex += fragment[2:]
                    sta = dex[0] if dex[0] > fragment[0] else fragment[0]
                    end = dex[1] if dex[1] < fragment[1] else fragment[1]
                    new_fragment = [sta, end] + fragment[2:] + dex_info
                    mapped_fragment.append(new_fragment)
                    if fragment[1] > dex[1]:
                        tmp_fragment.append([dex[1],
                                            fragment[1]] + fragment[2:])
        else:
            if flag:
                return mapped_fragment
            else:
                return index

import re

def ismember(a, b):
    """ https://stackoverflow.com/questions/15864082/python-equivalent-of-matlabs-ismember-function?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
    """
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = i
    return [bind.get(itm, None) for itm in a]  # None can be replaced by any other "not in b" value
      
def alphanum_key(s):
    """
    A function that gets the value of integers included in the string s.
    From StackOverflow by jamylak (http://stackoverflow.com/users/1219006/jamylak).
    http://stackoverflow.com/questions/11339210/how-to-get-integer-values-from-a-string-in-python/11339230#11339230
    """
    out = map(int, re.findall(r'\d+', s))
    out = list(out)
    return out[0]    
    
def find_by_num(lst, node_nums):
    """
    A function that searches lst (1st argument, it should be a list/tuple of objects with  
    a loc_name property) for an object containing the numbers from node_nums (2nd argument, 
    it can be an int or a list or tuple of ints) in its loc_name property. 
    The function returns the object if it is found, else None is returned.
    """
    obj = None
    node_nums = set(node_nums)
    for item in lst:
        item_nums = alphanum_key(item.loc_name)        
        try:
            item_nums = set(item_nums)
            if item_nums == node_nums:
                obj = item
                break             
        except:
            continue
    return obj    

def strsci2strfix(contents):
    """
    This should replace the occurencies of numeric strings in scientific notation with
    strings corresponding to the same numbers, but expressed in fixed-point notation.
    The underlying regexp to filter numbers in scientific notation is from
    https://stackoverflow.com/questions/18152597/extract-scientific-number-from-string
    and it should hopefully work well in the majority of situations.
    """
    sci_num = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?')
    for idx,_ in enumerate(contents):
        line = contents[idx]
        hexnum = line.find('#')>=0 or line.find('"')>=0
        if hexnum:        
            continue
        strsci = [x for x in re.findall(sci_num,line)]        
        for scistr in strsci:
            scistr = scistr.replace(' ','')
            scistr = str(float(scistr))
            fixstr = ('%f' % float(scistr))
            whole, fract = fixstr.split('.')
            if float(fract)>0:
                fixstr = fixstr.rstrip('0').rstrip('.')
            else:
                fixstr = whole+'.0'
            line = line.replace(scistr,fixstr)
        contents[idx] = line   
    return contents    

def recursive_check(fcn, inp, ref, max_depth):
    """
    The function to apply a custom function fcn (the 1st input) on each element of inp.
       
    fcn should take two inputs: inp (the 2nd input) and ref (the 3rd input). 
    inp can be a scalar or a (nested with many nested levels) list or tuple or a range.  
    
    The fcn should "compare" each entry of inp with ref and return True or False. 
    The the overall output of recursive_check is a conjunction of these True/False answers 
    for each checked element.
    
    The maximal number of nested levels to consider can be modified by passing the 4th 
    input argument (a nonnegative integer).
    """
    if max_depth==0: 
        return fcn(inp, ref)
    elif isinstance(inp, (list, tuple, range)):
        go_down = [recursive_check(fcn, item, ref, max_depth-1) for item in inp]
        return all(go_down)
    else:
        return fcn(inp, ref)


def flatten(inp):
    """
    The function to flatten a nested list or tuple of lists, tuples or ranges,
    with any number of nested levels.
    From StackOverflow by James Brady (http://stackoverflow.com/users/29903/james-brady)
    http://stackoverflow.com/questions/406121/flattening-a-shallow-list-in-python/406822#406822
    """
    out = []
    for el in inp:
        if isinstance(el, (list, tuple, range)):
            out.extend(flatten(el))
        else:
            out.append(el)
    return out
    
    
def check_contents(inp, typ, max_depth = 99):
    """
    The function to check inp (the 1st input) to be of a certain class, or check all its 
    contents (including contents of nested lists, tuples, ranges) against a certain class 
    in case inp is a list, tuple, range.
    The 2nd input provides the class to be checked, and its own class should be type.
    The maximal number of nested levels to consider can be modified by passing the 3rd 
    input argument (a positive integer).
    """    
    if not isinstance(max_depth, int) or max_depth < 0:
        raise TypeError('[pyUtils.check_contents] The 3rd input should be a positive int')
    fcn = isinstance
    out = recursive_check(fcn, inp, typ, max_depth)
    return out
        
        
def check_comparison(inp, op, cut = 0):
    """
    A function to compare inp (the 1st input) with a certain threshold, or compare all  
    its contents (including contents of nested lists, tuples, ranges) against a certain 
    threshold in case inp is a list, tuple, range.
    The 2nd input is the comparison operation, it should be a string ('>', '<', '>=',  
    '<=' or '=='). 
    The 3nd input provides the threshold to be checked, and its own class should be int
    or float. The default threshold is zero (threrefore the 3rd input is optional).    
    """
    err_3rd_inp = ('[pyUtils.check_comparison] The 3rd input should be either int or '
                   'float, or omitted (then the default value is zero)') 
    if not isinstance(cut, (int, float)):
        raise TypeError(err_3rd_inp)
    cond = check_contents(inp, (int, float))  # first check if comparison makes sense
    if not cond:
        raise TypeError('[pyUtils.check_comparison] The 1st input should be (real) numeric')
    # Implement comparison
    import operator
    ops = {'>': operator.gt,
           '<': operator.lt,
           '>=': operator.ge,
           '<=': operator.le,
           '==': operator.eq}    
    fcn = ops[op]
    inp = flatten(inp)  # make a flat list
    return all(fcn(item, cut) for item in inp)


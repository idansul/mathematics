from sympy import symbols
import re
import matplotlib.pyplot as plt

x, D, a, b = symbols("x, D, a, b")


def elements_search_histogram(algebra, left_removals=[None], right_removals=[None]):
    # Finds the single brackets of all the elements in the algebra and returns them with a histogram of the steps and the increase indexes
    histogram = {0:0}
    elements = []
    increase_indexes = []
    cardinality = algebra.p ** algebra.dimension
    for i in range(1, cardinality):
        chunk = algebra.get_all_elements_dict(start=i, stop=i+1, left_removals=left_removals, right_removals=right_removals)
        chunk = list(set(chunk.values()) - set(elements))
        elements.extend(chunk)
        histogram[i] = len(elements)
        if histogram[i] - histogram[i-1]:
            increase_indexes.append(i)
        print(i, len(elements))
        if len(elements) == cardinality:
            break
    return elements, histogram, increase_indexes

def search_skip(i, p, hist):
    # Skips every p elements or when new elements are found in the search process
    if i == 1:
        return False
    elif hist[i] - hist[i-1] != 0:
        return False
    elif i % p != 0:
        return True
    return False

def esh_plot(esh, plot_type="plot"):
    # Plots the histogram of the elements search
    index, count = zip(*esh.items())
    eval(f"plt.{plot_type}(index, count)")
    plt.title("Brackets search histogram")
    plt.xlabel("Iterations")
    plt.ylabel("Number of elements")
    plt.show()
    
def is_component(components, element):
    # Checks if the element is composed of one of the components
    if element == 0:
        return False
    element = list(map(eval, re.split('[+-]', str(element))))
    element = [re.findall("D.*", str(value))[0] for value in element]
    element = list(map(eval, element))
    components = [comp * D for comp in components]
    if set(components) <= set(element):
        return True
    return False

def elements_with_components_of_length_l(elements, components, l=0):
    # Finds the elements with the given components and with length l and counts them
    count = 0
    for item in elements.items():
        if is_component([x, x**4], item[1]):
            if not l or l and len(re.split("[+]", str(item[1]))) == l:
                print(item)
                count += 1
    print(f"Number of elements: {count}")

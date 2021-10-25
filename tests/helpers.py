from typing import Any


def assert_objects(left: Any, right: Any, attr: str) -> bool:
    """
    Asserts if two objects are from the same class and have equal `attrs`.

    Parameters:
        left (Any): an object to compare with `right`.
        right (Any): an object to compare with `left`.
        attr (str): a string representing the attribute to compare between the objects.

    Returns:
        bool: `True` if objects are the same type AND have the same value for `attr`,
              else `False`.
    """
    assert hasattr(left, attr), f"Object `left` has no attribute {attr}."
    assert hasattr(right, attr), f"Object `right` has no attribute {attr}."

    if type(left) is type(right):
        return getattr(left, attr) == getattr(right, attr)

    return False


def batch_assert_objects(
    left: list[Any], right: list[Any], attr: str, assume_order: bool = False
) -> bool:
    """
    Iterates over two lists of objects and compare if both lists are the same.
    Calls `assert_objects` internally, so it evaluates if objects are the same type and
    have the same value for attr. By default, assumes lists aren't in the same order and
    orders them first.

    Parameters:
        left (list[Any]): a list of objects to compare with `right`.
        right (list[Any]): a list of objects to compare with `left`.
        attr (str): a string representing the attribute to compare between the objects.
        assume_order (bool): wheter to assume lists are ordered or not. If `False`,
                             orders the lists first. Default: False.

    Returns:
        bool: `True` if all objects in both lists, in order, are the same type and have
              the same value for `attr`, else `False`.
    """
    assert len(left) == len(right), "Lists aren't the same length."

    if not assume_order:
        left.sort()
        right.sort()

    for left_item, right_item in zip(left, right):
        if assert_objects(left_item, right_item, attr) == False:
            return False

    return True

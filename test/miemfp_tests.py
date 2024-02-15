from typing import List, Tuple
import asyncio
import numpy as np
from hypothesis import errors, given, settings, strategies as st # type: ignore

def compare_bhmie_functions(f1, f2, event1: asyncio.Event, event2: asyncio.Event) -> Tuple[bool, str]:
    async def async_closure(
        x: float,
        cxref: Tuple[float, float],
        cxs1: List[Tuple[float, float]],
        cxs2: List[Tuple[float, float]],
    ) -> bool:
        cxref = complex(*cxref)
        cxs1 = [complex(*c) for c in cxs1]
        cxs2 = [complex(*c) for c in cxs2]
        
        # This is to ensure that only one instance of each function is running at a time
        # to avoid memory issues in the FFI code
        await event1.wait()
        f1_result = f1(x, cxref, 2, cxs1, cxs2)[:2]
        await event2.wait()
        f2_result = f2(x, cxref, 2, cxs1, cxs2)[:2]
        
        return np.all(np.isclose(f1_result, f2_result))
    
    @settings(deadline=None)
    @given(
        # Must be bigger than an atom but still nanoscale
        x=st.floats(min_value=0.1, max_value=100),
        # Refractive indeces must be within a physically reasonable range
        cxref=st.tuples(st.floats(min_value=0.1, max_value=4.0), st.floats(min_value=0.1, max_value=4.0)),
        cxs1=st.lists(st.tuples(st.floats(min_value=0.1, allow_infinity=False), st.floats(min_value=0.1, allow_infinity=False)), min_size=10, max_size=100),
        cxs2=st.lists(st.tuples(st.floats(min_value=0.1, allow_infinity=False), st.floats(min_value=0.1, allow_infinity=False)), min_size=10, max_size=100),
    )
    def sync_closure(
        x: float,
        cxref: Tuple[float, float],
        cxs1: List[Tuple[float, float]],
        cxs2: List[Tuple[float, float]],
    ) -> bool:
        assert asyncio.run(async_closure(x, cxref, cxs1, cxs2))
    
    try:
        sync_closure()
        return True, "Test passed"
    except AssertionError as e:
        return False, f"AssertionError: {str(e)}"
    except errors.HypothesisException as e:
        return False, f"HypothesisException: {str(e)}"
import functools
import itertools
import dataclasses
import fractions
import math
import typing
import json

import pint

Numeric = typing.Union[float, int, fractions.Fraction]

ureg = pint.UnitRegistry()


@functools.total_ordering
@dataclasses.dataclass()
class Thread:
    """A thread

    Args:
        major_diameter: The thread's major diameter
        pitch: The thread's pitch; if not provided, will be calculated from the
            tpi
        tpi: The thread's TPI; if not provided, will be calculated from
            the pitch
    """

    major_diameter: pint.Quantity
    minor_diameter: pint.Quantity = dataclasses.field(init=False)

    pitch: pint.Quantity = None
    tpi: pint.Quantity = None

    def __post_init__(self) -> None:
        if self.pitch:
            self.tpi = 1 / self.pitch
        elif self.tpi:
            self.pitch = 1 / self.tpi
        else:
            raise ValueError

        self.minor_diameter = self.major_diameter - (
            (5 * math.sqrt(3) * self.pitch) / 8
        )

    def __gt__(self, other: "Thread") -> bool:
        """Order threads based first on major diameter & then on pitch"""
        return (self.major_diameter, self.minor_diameter) > (
            other.major_diameter,
            other.minor_diameter,
        )


class ISOThread(Thread):
    """A Thread following the ISO standard"""

    def __init__(self, major_diameter: pint.Quantity, pitch: pint.Quantity) -> None:
        """
        Args:
            major_diameter: The thread's major diameter
            pitch: The thread's pitch
        """
        super().__init__(major_diameter=major_diameter, pitch=pitch)

    @classmethod
    def from_mm(cls, major_diameter: Numeric, pitch: Numeric) -> "ISOThread":
        """Construct an ISO thread, given values in milimeters.

        Args:
            major_diameter: The thread's major diameter, in milimeters
            pitch: The thread's pitch, in milimeters
        """
        return cls(major_diameter * ureg.mm, pitch * ureg.mm)

    def __str__(self) -> str:
        return f"M{self.major_diameter.to(ureg.mm).magnitude}-{self.pitch.to(ureg.mm).magnitude}"


class UTSThread(Thread):
    """A thread follpowing the Unified Thread Standard"""

    def __init__(self, major_diameter: pint.Quantity, tpi: pint.Quantity) -> None:
        """
        Args:
            major_diameter: The thread's major diameter
            tpi: The thread's teeth per inch
        """
        super().__init__(major_diameter=major_diameter, tpi=tpi)
        self.formatted_tpi = f"{self.tpi.to(1 / ureg.inch).magnitude}"


class FractionalUTSThread(UTSThread):
    """A UTS thread specified by a fractional major diameter"""

    @classmethod
    def from_fraction(
        cls, numerator: int, denominator: int, tpi: int
    ) -> "FractionalUTSThread":
        """Construct a UTS Thread, given values in fractional inches

        Args:
            numerator: The numerator of the thread's major diameter, in inches
            denominator: The denominator of the thread's major diameter, in inches
            tpi: The thread's teeth per inch
        """
        return cls(
            fractions.Fraction(numerator, denominator) * ureg.inch, tpi / ureg.inch
        )

    def __str__(self) -> str:
        return f'{self.major_diameter.to(ureg.inch).magnitude}"-{self.formatted_tpi}'


class NumberedUTSThread(UTSThread):
    """A UTS thread specified by a numbered major diameter"""

    def __init__(
        self, major_diameter: pint.Quantity, tpi: pint.Quantity, number: int
    ) -> None:
        """
        Args:
            major_diameter: The thread's major diameter
            tpi: The thread's teeth per inch
            number: The thread's number
        """
        self.number = number
        super().__init__(major_diameter=major_diameter, tpi=tpi)

    @classmethod
    def from_number(cls, number: int, tpi: int) -> "NumberedUTSThread":
        """Construct a Numbered UTS Thread, given values in inches

        Args:
            number: The thread's number
            tpi: The thread's teeth per inch
        """
        return cls(
            (number * 0.013 * ureg.inch) + (0.06 * ureg.inch), tpi / ureg.inch, number
        )

    def __str__(self) -> str:
        return f"#{self.number}-{self.formatted_tpi}"


@dataclasses.dataclass(order=True)
class DifferentialThread:
    """A differential thread pair

    Args:
        thread_pair: The pair of threads that define the differential
    """

    thread_pair: typing.Tuple[Thread, Thread] = dataclasses.field(compare=False)

    effective_pitch: pint.Quantity = dataclasses.field(init=False)
    effective_tpi: pint.Quantity = dataclasses.field(init=False, compare=False)

    radial_clearance: pint.Quantity = dataclasses.field(init=False)

    def __post_init__(self) -> None:
        smaller_thread, larger_thread = sorted(self.thread_pair)

        self.effective_pitch = abs(smaller_thread.pitch - larger_thread.pitch)
        self.effective_tpi = 1 / self.effective_pitch

        radial_clearance = (
            larger_thread.minor_diameter - smaller_thread.major_diameter
        ) / 2
        if radial_clearance.magnitude > 0:
            self.radial_clearance = radial_clearance
        else:
            self.radial_clearance = 0.0 * ureg.inch


THREADS = [
    ISOThread.from_mm(1, 0.25),
    ISOThread.from_mm(1, 0.2),
    ISOThread.from_mm(1.2, 0.25),
    ISOThread.from_mm(1.2, 0.2),
    ISOThread.from_mm(1.4, 0.3),
    ISOThread.from_mm(1.4, 0.2),
    ISOThread.from_mm(1.6, 0.35),
    ISOThread.from_mm(1.6, 0.2),
    ISOThread.from_mm(1.8, 0.35),
    ISOThread.from_mm(1.8, 0.2),
    ISOThread.from_mm(2, 0.4),
    ISOThread.from_mm(2, 0.25),
    ISOThread.from_mm(2.5, 0.45),
    ISOThread.from_mm(2.5, 0.35),
    ISOThread.from_mm(3, 0.5),
    ISOThread.from_mm(3, 0.35),
    ISOThread.from_mm(3.5, 0.6),
    ISOThread.from_mm(3.5, 0.35),
    ISOThread.from_mm(4, 0.7),
    ISOThread.from_mm(4, 0.5),
    ISOThread.from_mm(5, 0.8),
    ISOThread.from_mm(5, 0.5),
    ISOThread.from_mm(5.5, 0.9),
    ISOThread.from_mm(5.5, 0.5),
    ISOThread.from_mm(6, 1.0),
    ISOThread.from_mm(6, 0.75),
    ISOThread.from_mm(7, 1),
    ISOThread.from_mm(7, 0.75),
    ISOThread.from_mm(8, 1.25),
    ISOThread.from_mm(8, 1.0),
    ISOThread.from_mm(8, 0.75),
    ISOThread.from_mm(10, 1.5),
    ISOThread.from_mm(10, 1.25),
    ISOThread.from_mm(10, 1.0),
    ISOThread.from_mm(12, 1.75),
    ISOThread.from_mm(12, 1.5),
    ISOThread.from_mm(12, 1.25),
    ISOThread.from_mm(14, 2.0),
    ISOThread.from_mm(14, 1.5),
    ISOThread.from_mm(16, 2.0),
    ISOThread.from_mm(16, 1.5),
    ISOThread.from_mm(18, 2.5),
    ISOThread.from_mm(18, 2.0),
    ISOThread.from_mm(18, 1.5),
    ISOThread.from_mm(20, 2.5),
    ISOThread.from_mm(20, 2.0),
    ISOThread.from_mm(20, 1.5),
    ISOThread.from_mm(22, 2.5),
    ISOThread.from_mm(22, 2.0),
    ISOThread.from_mm(22, 1.5),
    ISOThread.from_mm(24, 3.0),
    ISOThread.from_mm(24, 2.0),
    ISOThread.from_mm(27, 3.0),
    ISOThread.from_mm(27, 2.0),
    ISOThread.from_mm(30, 3.5),
    ISOThread.from_mm(30, 2.0),
    ISOThread.from_mm(33, 3.5),
    ISOThread.from_mm(33, 2.0),
    ISOThread.from_mm(36, 4.0),
    ISOThread.from_mm(36, 3.0),
    ISOThread.from_mm(39, 4.0),
    ISOThread.from_mm(39, 3.0),
    ISOThread.from_mm(42, 4.5),
    ISOThread.from_mm(42, 3.0),
    ISOThread.from_mm(45, 4.5),
    ISOThread.from_mm(45, 3.0),
    ISOThread.from_mm(48, 5.0),
    ISOThread.from_mm(48, 3.0),
    ISOThread.from_mm(52, 5.0),
    ISOThread.from_mm(52, 4.0),
    ISOThread.from_mm(56, 5.5),
    ISOThread.from_mm(56, 4.0),
    ISOThread.from_mm(60, 5.5),
    ISOThread.from_mm(60, 4.0),
    ISOThread.from_mm(62, 6.0),
    ISOThread.from_mm(62, 4.0),
    NumberedUTSThread.from_number(0, 80),
    NumberedUTSThread.from_number(1, 64),
    NumberedUTSThread.from_number(1, 72),
    NumberedUTSThread.from_number(2, 56),
    NumberedUTSThread.from_number(2, 64),
    NumberedUTSThread.from_number(3, 48),
    NumberedUTSThread.from_number(3, 56),
    NumberedUTSThread.from_number(4, 40),
    NumberedUTSThread.from_number(4, 48),
    NumberedUTSThread.from_number(5, 40),
    NumberedUTSThread.from_number(5, 44),
    NumberedUTSThread.from_number(6, 32),
    NumberedUTSThread.from_number(6, 40),
    NumberedUTSThread.from_number(8, 32),
    NumberedUTSThread.from_number(8, 36),
    NumberedUTSThread.from_number(10, 24),
    NumberedUTSThread.from_number(10, 28),
    NumberedUTSThread.from_number(12, 24),
    NumberedUTSThread.from_number(12, 38),
    NumberedUTSThread.from_number(12, 32),
    FractionalUTSThread.from_fraction(1, 4, 20),
    FractionalUTSThread.from_fraction(1, 4, 28),
    FractionalUTSThread.from_fraction(1, 4, 32),
    FractionalUTSThread.from_fraction(5, 16, 18),
    FractionalUTSThread.from_fraction(5, 16, 24),
    FractionalUTSThread.from_fraction(5, 16, 32),
    FractionalUTSThread.from_fraction(3, 8, 16),
    FractionalUTSThread.from_fraction(3, 8, 24),
    FractionalUTSThread.from_fraction(3, 8, 32),
    FractionalUTSThread.from_fraction(7, 16, 14),
    FractionalUTSThread.from_fraction(7, 16, 20),
    FractionalUTSThread.from_fraction(7, 16, 28),
    FractionalUTSThread.from_fraction(1, 2, 13),
    FractionalUTSThread.from_fraction(1, 2, 20),
    FractionalUTSThread.from_fraction(1, 2, 28),
    FractionalUTSThread.from_fraction(9, 16, 12),
    FractionalUTSThread.from_fraction(9, 16, 18),
    FractionalUTSThread.from_fraction(9, 16, 24),
    FractionalUTSThread.from_fraction(5, 8, 11),
    FractionalUTSThread.from_fraction(5, 8, 18),
    FractionalUTSThread.from_fraction(5, 8, 24),
    FractionalUTSThread.from_fraction(3, 4, 10),
    FractionalUTSThread.from_fraction(3, 4, 16),
    FractionalUTSThread.from_fraction(3, 4, 20),
    FractionalUTSThread.from_fraction(7, 8, 9),
    FractionalUTSThread.from_fraction(7, 8, 14),
    FractionalUTSThread.from_fraction(7, 8, 20),
]

with open("differential_threads.json", "w", encoding="utf-8") as json_file_object:
    json.dump(
        [
            {
                "threads": [str(thread) for thread in differential.thread_pair],
                "radial_clearance": f"{differential.radial_clearance.to(ureg.inch):0.3f~}",
                "effective_pitch": f"{differential.effective_pitch.to(ureg.inch):0.5f~}",
                "effective_tpi": f"{differential.effective_tpi.to(1 / ureg.inch).magnitude:0.2f}",
            }
            for differential in sorted(
                (
                    DifferentialThread(thread_pair)
                    for thread_pair in (
                        (left_thread, right_thread)
                        for left_thread, right_thread in itertools.combinations(
                            THREADS, 2
                        )
                        if left_thread.pitch != right_thread.pitch
                    )
                )
            )
        ],
        json_file_object,
        ensure_ascii=False,
        indent=2,
    )

from dataclasses import dataclass

@dataclass
class Bench():
    thread_id: int
    task: str
    start: int
    end: int
    duration: int
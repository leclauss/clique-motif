from multiprocessing import Process, Queue


class ParallelRun:
    def __init__(self, function, threads, outputFunction=None, log=lambda *_, **__: None):
        self.function = function
        self.log = log
        self.outputFunction = outputFunction

        self.argsQueue = Queue()
        self.workers = [Process(target=self.__work, args=(i,)) for i in range(threads)]
        if outputFunction is not None:
            self.resultQueue = Queue()
            self.outputWorker = Process(target=self.__outputWorker)
            self.outputWorker.start()
        for worker in self.workers:
            worker.start()

    def __work(self, threadId):
        while True:
            args = self.argsQueue.get()
            if args is None:
                self.argsQueue.put(None)
                break
            self.log("Thread", threadId, "running", *args)
            result = self.function(*args)
            if self.outputFunction is not None:
                self.resultQueue.put(result)
        self.log("Thread", threadId, "done")

    def __outputWorker(self):
        while True:
            result = self.resultQueue.get()
            if result is None:
                break
            self.outputFunction(*result)

    def run(self, args):
        self.argsQueue.put(args)

    def join(self):
        self.argsQueue.put(None)
        for worker in self.workers:
            worker.join()
        if self.outputFunction is not None:
            self.resultQueue.put(None)
            self.outputWorker.join()

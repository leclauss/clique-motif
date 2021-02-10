from multiprocessing import Process, Queue


def parallelRun(function, argsQueue, threads):
    def __work(func, inputQueue, outputQueue):
        while True:
            try:
                args = inputQueue.get_nowait()
            except Exception:
                break
            outputQueue.put(func(*args))

    resultQueue = Queue()
    workers = [Process(target=__work, args=(function, argsQueue, resultQueue)) for _ in range(threads)]
    for worker in workers:
        worker.start()
    for worker in workers:
        worker.join()
    return resultQueue

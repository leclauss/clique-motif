# TODO

header = "Time Series Length,Method,Motif Shape,Noise,Motif Size,Window,Actual Range,Input Range," \
         "Algorithm1,Runtime,Precision,Recall,F1,Found Size,..."  # scores for generator motif if input = actual range


# also always scores for ground truth = cliqueMotif

# motif = [int(index) for index in info["matchings"].split(",")]


def getPointBasedScores(motif, motifsFound, window):
    # get relevant motif points
    motifPoints = getMotifPoints(motif, window)

    bestStats = (0.0, 0.0, 0.0, 0, -1)  # f1Score, precision, recall, foundSize, index
    # calculate measures for each motif found
    for index, motifFound in enumerate(motifsFound):
        foundPoints = getMotifPoints(motifFound, window)
        truePositives = len(motifPoints.intersection(foundPoints))
        falsePositives = len(foundPoints.difference(motifPoints))
        falseNegatives = len(motifPoints.difference(foundPoints))
        # calculate precision
        returnedPositives = truePositives + falsePositives
        precision = (truePositives / returnedPositives if returnedPositives > 0 else 0.0)
        # calculate recall
        actualPositives = truePositives + falseNegatives
        recall = (truePositives / actualPositives if actualPositives > 0 else 0.0)
        # calculate f1 score
        f1Score = (2 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0)
        # save best f1 score
        if f1Score > bestStats[0]:
            bestStats = (f1Score, precision, recall, len(motifFound), index)
    return bestStats


def getMotifPoints(motif, window):
    motifPoints = set()
    for index in motif:
        motifPoints = motifPoints.union(set(range(index, index + window)))
    return motifPoints

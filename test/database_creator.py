import os
import sqlite3
import pandas as pd
import numpy as np
from sqlalchemy import Column, String, Integer, REAL, create_engine
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import sessionmaker

# %%
Base = declarative_base()


class PairwisePrediction(Base):
    __tablename__ = 'AF2 full-lib'

    pid = Column("pid", String, primary_key=True)
    data_path = Column("data_path", String)
    geneA = Column("gene A", String)
    geneB = Column("gene B", String)
    nameA = Column("name A", String)
    nameB = Column("name B", String)
    seqA = Column("seq A", String)
    seqB = Column("seq B", String)
    canonicA = Column("canonic A", Integer)
    canonicB = Column("canonic B", Integer)
    lenA = Column("len A", Integer)
    lenB = Column("len B", Integer)
    pTM = Column("pTM", REAL)
    ipTM = Column("ipTM", REAL)
    pDockQ = Column("pDockQ", REAL)
    meanpDockQ = Column("meanpDockQ", REAL)

    def __init__(self,
                 pid,
                 data_path,
                 geneA,
                 geneB,
                 nameA,
                 nameB,
                 seqA,
                 seqB,
                 canonicA,
                 canonicB,
                 lenA,
                 lenB,
                 pTM,
                 ipTM,
                 pDockQ,
                 meanpDockQ, ):
        self.pid = pid
        self.meanpDockQ = meanpDockQ
        self.pDockQ = pDockQ
        self.ipTM = ipTM
        self.pTM = pTM
        self.lenB = lenB
        self.lenA = lenA
        self.canonicB = canonicB
        self.canonicA = canonicA
        self.seqB = seqB
        self.seqA = seqA
        self.nameB = nameB
        self.nameA = nameA
        self.geneB = geneB
        self.data_path = data_path
        self.geneA = geneA

    def __repr__(self):
        return f"{self.__tablename__}({self.pid}, {self.data_path}, {self.geneA}, {self.geneB}, {self.nameA}, {self.nameB}, {self.seqA}, {self.seqB}, {self.canonicA}, {self.canonicB}, {self.lenA}, {self.lenB}, {self.pTM}, {self.ipTM}, {self.pDockQ}, {self.meanpDockQ})"


# %%
engine = create_engine('sqlite:///PPIDB_full.db', echo=True)
Base.metadata.create_all(bind=engine)

Session = sessionmaker(bind=engine)
session = Session()

interaction = PairwisePrediction(
    pid=str(hex(np.random.randint(10 ** 10))),
    data_path='path/to/data',
    geneA='geneA',
    geneB='geneB',
    nameA='nameA',
    nameB='nameB',
    seqA='seqA',
    seqB='seqB',
    canonicA=1,
    canonicB=1,
    lenA=100,
    lenB=100,
    pTM=0.5,
    ipTM=0.5,
    pDockQ=0.5,
    meanpDockQ=0.5)

print(interaction)

# %%
# session.add(interaction)
# session.commit()


fig


from peewee import *

db = SqliteDatabase('RMF.db')

class RMF(Model):
    class Meta:
        database = db

class RMFClasses(RMF):
    class__ = CharField()

class RMFParams(RMFClasses):
    name = CharField()
    type = ForeignKeyField(RMFClasses, related_name='params')
    module = CharField()
    constants = CharField()
    massfile = CharField()
    tabfile = CharField()
    
db.connect()
RMFParams.create_table(fail_silently=True)
RMFClasses.create_table(fail_silently=True)


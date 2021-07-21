"""
CBMPy: CBNetDB module
=====================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2022 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

Author: Brett G. Olivier
Contact email: bgoli@users.sourceforge.net
Last edit: $Author: bgoli $ ($Id: CBNetDB.py 710 2020-04-27 14:22:34Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import

# from __future__ import unicode_literals

import os
import time
import re
import webbrowser
import csv

try:
    from urllib2 import quote as urlquote
    from urllib2 import unquote as urlunquote
except ImportError:
    from urllib.parse import quote as urlquote
    from urllib.parse import unquote as urlunquote

HAVE_SQLITE2 = False
HAVE_SQLITE3 = False

try:
    import sqlite3 as sqlite3

    HAVE_SQLITE3 = True
except ImportError:
    HAVE_SQLITE3 = False
    try:
        from pysqlite2 import dbapi2 as sqlite2

        HAVE_SQLITE2 = True
    except ImportError:
        HAVE_SQLITE2 = False

from .CBConfig import __CBCONFIG__ as __CBCONFIG__

__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']


class NetDBbase(object):
    text_encoding = 'utf8'

    def URLEncode(self, txt):
        """
        URL encodes a string.

        """
        return urlquote(txt.encode(self.text_encoding))

    def URLDecode(self, txt):
        """
        Decodes a URL encoded string

        """
        return urlunquote(txt)


class DBTools(NetDBbase):
    """
    Tools to work with SQLite DB's (optimized, no SQL required).

    """

    sqlite = None
    sqlite_version = None
    db_conn = None
    db_cursor = None
    db_tables = None

    def __init__(self):
        if HAVE_SQLITE3:
            self.sqlite = sqlite3
            self.sqlite_version = 3
        elif HAVE_SQLITE2:
            self.sqlite = sqlite2
            self.sqlite_version = 2
        else:
            raise RuntimeError("\nSQLite not installed")

        self.db_tables = []

    def connectSQLiteDB(self, db_name, work_dir=None):
        """
        Connect to a sqlite database.

        - *db_name* the name of the sqlite database
        - *work_dir* the optional database path

        """
        # connect to DB
        if work_dir != None:
            self.db_conn = self.sqlite.connect(os.path.join(work_dir, db_name))
        else:
            self.db_conn = self.sqlite.connect(db_name)

        self.db_cursor = self.db_conn.cursor()

    def commitDB(self):
        """
        Commits all curent changes to DB, returns a boolean.
        """
        try:
            self.db_cursor.connection.commit()
            return True
        except:
            return False

    def createDBTable(self, table, sqlcols):
        """
        Create a database table if it does not exist:

         - *table* the table name
         - *sqlcols* a list containing the SQL definitions of the table columns: <id> <type> for example `['gene TEXT PRIMARY KEY', 'aa_seq TEXT', 'nuc_seq TEXT', 'aa_len INT', 'nuc_len INT']`

        Effectively writes CREATE TABLE "table" (<id> <type>, gene TEXT PRIMARY KEY, aa_seq TEXT, nuc_seq TEXT, aa_len INT, nuc_len INT) % table
        """
        SQL = 'CREATE TABLE %s (' % table
        for c in sqlcols:
            SQL += ' %s,' % c
        SQL = SQL[:-1]
        SQL += ' )'
        print(SQL)
        try:
            self.db_cursor.execute('SELECT * FROM %s' % table)
            print('Table {} exists'.format(table))
        except self.sqlite.OperationalError:
            print('Table {} does not exist, creating it'.format(table))
            self.db_cursor.execute(SQL)
        self.db_tables.append(table)

    def insertData(self, table, data, commit=True):
        """
        Insert data into a table: "INSERT INTO %s (?, ?, ?, ?, ?) VALUES (?, ?, ?, ?, ?)" % table,
                                    (?, ?, ?, ?, ?)) )

         - *table* the DB table name
         - *data* a dictionary of {id:value} pairs
         - *commit* whether to commit the data insertions

        """
        colstr = "("
        valstr = "VALUES ("
        vals = []
        # I want to use dictionaries but this just keeps things backwards compatible for a while
        try:
            data.keys()
            for d in data:
                colstr += '%s, ' % d
                valstr += '?, '
                vals.append(data[d])
        except AttributeError:
            print(
                '\n\nWARNING: data now uses a dictionary as input please update your code - see docstring for details\n'
            )
            time.sleep(3)
            for d in data:
                colstr += '%s, ' % d[0]
                valstr += '?, '
                vals.append(d[1])
        colstr = colstr[:-2] + ')'
        valstr = valstr[:-2] + ')'
        sql = "INSERT INTO %s %s %s" % (table, colstr, valstr)
        # print sql
        # time.sleep(1)
        try:
            self.db_cursor.execute(sql, vals)
            if commit:
                self.db_cursor.connection.commit()
            return True
        except AttributeError:
            return False

    def updateData(self, table, col, rid, data, commit=True):
        """
        Update already defined data

         - *table* the table name
         - *col* the column name
         - *rid* the row id to update
         - *data* a dictionary of {id:value} pairs
         - *commit* whether to commit the data updates

         UPDATE COMPANY SET ADDRESS = 'Texas' WHERE ID = 6;

        """

        sql = 'UPDATE {} SET '.format(table)
        for d in data:
            sql += '{}=\"{}\", '.format(d, str(data[d]).replace('\"', '\''))
        sql = sql[:-2] + ' WHERE {}=\"{}\"'.format(col, rid)

        # print(sql)
        try:
            self.db_cursor.execute(sql)
            if commit:
                self.db_cursor.connection.commit()
            return True
        except AttributeError:
            return False

    def checkEntryInColumn(self, table, col, rid):
        """
        Check if an entry exists in a table

        - *table* the table name
        - *col* the column name
        - *rid* the row to search for

        """
        self.db_cursor.execute(
            "SELECT count(*) FROM {} WHERE {}=\"{}\"".format(table, col, rid)
        )
        data = self.db_cursor.fetchone()[0]
        if data == 0:
            # print('There is no component named {}'.format(rid))
            return False
        else:
            # print('Component {} found in {} row(s)'.format(rid, data))
            return True

    def executeSQL(self, sql):
        """
        Execute a SQL command:

         - *sql* a string containing a SQL command

        """
        try:
            self.db_cursor.execute(sql)
            return True
        except AttributeError as ex:
            print('Error executing command')
            print(ex)
            return False

    def getColumns(self, table, cols):
        """
        Fetch the contents of one or more columns of data in a table

         - *table* the database table
         - *cols* a list of one or more column id's

        """
        sql = "SELECT "
        for c in cols:
            sql += '{}, '.format(str(c))
        sql = sql[:-2]
        sql += ' FROM {}'.format(table)
        dout = [[] for i in range(len(cols))]
        # print(sql)
        # print(dout)
        try:
            data = self.db_cursor.execute(sql).fetchall()
            for r_ in data:
                for c_ in range(len(r_)):
                    dout[c_].append(str(r_[c_]))
            del data
        except AttributeError:
            return None
        else:
            return dout

    def getRow(self, table, col, rid):
        """
        Get the table row(s) which correspond to rid in column. Returns the row(s) as a list, if the column is the primary key
        this is always a single entry.

         - *table* the database table
         - *col* the column id
         - *rid* the row index id

        """
        sql = "SELECT * FROM {} WHERE {}=\"{}\"".format(table, col, rid)
        # print(sql)
        dout = []
        try:
            data = self.db_cursor.execute(sql).fetchall()
            if len(data) > 1:
                print(
                    'INFO: getRow is returning multiple rows for query id: {}'.format(
                        rid
                    )
                )
            for r_ in data:
                dout.append([str(d) for d in r_])
            del data
        except AttributeError:
            return None
        else:
            return dout

    def getCell(self, table, col, rid, cell):
        """
        Get the table cell which correspond to rid in column. Returns the value or None

         - *table* the database table
         - *col* the column id
         - *rid* the row index id
         - *cell* the column of the cell you want tp extract

        """
        sql = "SELECT {} FROM {} WHERE {}=\"{}\"".format(cell, table, col, rid)
        # print(sql)
        data = None
        try:
            data = self.db_cursor.execute(sql).fetchone()
            if data is None:
                return None
            data = str(self.db_cursor.execute(sql).fetchone()[0])
        except AttributeError:
            return None
        return data

    def getTable(self, table, colOut=False):
        """
        Returns an entire database table

         - *table* the table name
         - *colOut* optionally return a tuple of (data,ColNames)

        """
        sql = 'SELECT * FROM %s' % table
        sql2 = "PRAGMA table_info( %s )" % table

        r = None
        col = None
        try:
            r = self.db_cursor.execute(sql).fetchall()
            if colOut:
                col = self.db_cursor.execute(sql2).fetchall()
                col = [str(a[1]) for a in col]
        except (AttributeError, sqlite3.OperationalError):
            return None
        if colOut:
            return r, col
        else:
            return r

    def dumpTableToTxt(self, table, filename):
        """
        Save a table as tab separated txt file

         - *table* the table to export
         - *filename* the filename of the table dump

        """
        data, head = self.getTable(table, colOut=True)
        data.insert(0, head)
        from .CBTools import exportLabelledLinkedList

        exportLabelledLinkedList(data, fname=filename, names=None, sep='\t')

    def dumpTableToCSV(self, table, filename):
        """
        Save a table as tab separated txt file

         - *table* the table to export
         - *filename* the filename of the table dump

        """
        data, head = self.getTable(table, colOut=True)
        data.insert(0, head)

        F = open(filename, 'w')
        csvw = csv.writer(F, dialect='excel')
        csvw.writerows(data)
        F.close()
        del csvw

    def fetchAll(self, sql):
        """Raw SQL query e.g. 'SELECT id FROM gene WHERE gene=\"G\"' """
        # print(sql)
        r = None
        try:
            r = self.db_cursor.execute(sql).fetchall()
        except AttributeError as ex:
            print(ex)
        return r

    def closeDB(self):
        """Close the DB connection and reset the DBTools instance (can be reconnected)"""
        self.db_conn.close()
        self.db_conn = None
        self.db_cursor = None
        self.db_tables = []


class KeGGTools(object):
    """
    Class that holds useful methods for querying KeGG via a SUDS provided soap client
    """

    Kclient = None

    def __init__(self, url):
        import suds

        self.Kclient = suds.client.Client(url)

    def fetchSeqfromKeGG(self, k_gene):
        """
        Given a gene name try and retrieve the gene and amino acid sequence
        """
        g2 = 'None'
        p2 = 'None'
        try:
            g = self.Kclient.service.bget("-f -n n %s" % k_gene)
            if g == None:
                print(
                    '\n*****\nWARNING: potential naming error in gene: {}!!\n*****\n'.format(
                        k_gene
                    )
                )
            g2 = g.split('(N)')[1].replace('\n', '')
        except AttributeError as ex:
            print('\nGene sequence get exception ({})!\n'.format(k_gene))
            print(ex)
        try:
            p = self.Kclient.service.bget("-f -n a %s" % k_gene)
            if p == None:
                print(
                    '\n*****\nWARNING: potential naming error in gene: {}!!\n*****\n'.format(
                        k_gene
                    )
                )
            p2 = p.split('(A)')[1].replace('\n', '')
        except AttributeError as ex:
            print('\nProtein sequence get exception ({})!\n'.format(k_gene))
            print(ex)
        return g2, p2


class KeGGSequenceTools(object):
    """
    Using the KeGG connector this class provides tools to construct an organims specific sequence database
    """

    DB = None
    KEGG = None

    def __init__(self, url, db_name, work_dir):
        self.DB = DBTools()
        self.DB.connectSQLiteDB(db_name, work_dir)
        self.KEGG = KeGGTools(url)

    def buildGeneDatabase(
        self, genes, tablename, UPDATE_IF_EXISTS=False, default_length=0
    ):
        cntr = 1
        cntr2 = 1

        for ecg in genes:
            print('Processing gene {} of {}'.format(cntr, len(genes)))
            entry_exists = False
            testg = self.DB.db_cursor.execute(
                'SELECT * FROM %s WHERE gene="%s" ' % (tablename, ecg)
            ).fetchall()
            if len(testg) > 0:
                entry_exists = True

            # if ecg in ['eco:b1898','eco:b1899','eco:b3692','eco:b3111','eco:b4228','eco:b2978',\
            # 'eco:b1416','eco:b3112','eco:b1417','eco:b3768','eco:b3767','eco:b4229']:
            # raw_input(testg)

            tstart = time.time()
            if not entry_exists:
                print('\tadding gene {}'.format(ecg), end=" ")
                gene2, prot2 = self.KEGG.fetchSeqfromKeGG(ecg)
                if gene2 != 'None' and prot2 != 'None':
                    self.DB.db_cursor.execute(
                        "INSERT INTO %s (gene, aa_seq, nuc_seq, aa_len, nuc_len) VALUES (?, ?, ?, ?, ?)"
                        % tablename,
                        (
                            str(ecg),
                            str(prot2),
                            str(gene2),
                            int(len(prot2)),
                            int(len(gene2)),
                        ),
                    )
                else:
                    print(
                        '\nGene {} cannot be found and is probably an incorrect annotation assigning length: {}\n'.format(
                            ecg, default_length
                        )
                    )
                    self.DB.db_cursor.execute(
                        "INSERT INTO %s (gene, aa_seq, nuc_seq, aa_len, nuc_len) VALUES (?, ?, ?, ?, ?)"
                        % tablename,
                        (str(ecg), 'None', 'None', default_length, default_length),
                    )
            elif entry_exists and UPDATE_IF_EXISTS:
                print('\tupdating gene {}'.format(ecg), end=" ")
                gene2, prot2 = self.KEGG.fetchSeqfromKeGG(ecg)
                self.DB.db_cursor.execute(
                    'UPDATE %s SET aa_seq="%s", nuc_seq="%s", aa_len="%s", nuc_len="%s" WHERE gene="%s"'
                    % (tablename, prot2, gene2, int(len(prot2)), int(len(gene2)), ecg)
                )
            else:
                print('\tskipping gene {}'.format(ecg), end=" ")
            tend = time.time()
            print(' ... done ({}).'.format(tend - tstart))
            #  if cntr == 6:
            #  break
            cntr += 1
            cntr2 += 1
            if cntr2 == 21:
                self.DB.db_cursor.connection.commit()
                cntr2 = 1

    def getPeptideLengthsFromDB(self, genes, keg_prefix):
        gene_peplen = {}
        for G in genes:
            print(G)
            ##  Glen = self.cursor.execute('SELECT aa_len FROM gene_data WHERE gene="%s"' % G).fetchall()[0][0]
            Glen = self.DB.fetchAll('SELECT aa_len FROM gene_data WHERE gene="%s"' % G)[
                0
            ][0]
            print(Glen)
            gene_peplen.update({G.replace(keg_prefix, ''): Glen})
        return gene_peplen


class RESTClient(NetDBbase):
    """
    Class that provides the basis for application specific connectors to REST web services
    """

    site_root = None
    conn = None
    history = ''
    CONNECTED = False
    USER_AGENT = (
        'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:13.0) Gecko/20100101 Firefox/13.0'
    )

    def __init__(self):
        self.urllib2 = urllib2

    def Log(self, txt):
        """
        Add txt to logfile history

         - *txt* a string
        """
        self.history += '%s - %s\n' % (time.strftime('%H:%M:%S'), str(txt))

    def GetLog(self):
        """
        Return the logged history
        """
        return self.history

    def Connect(self, root):
        """
        Establish HTTP connection to

         - *root* the site root "www.google.com"

        """
        try:
            self.site_root = root
            self.conn = self.urllib2.httplib.HTTPConnection(self.site_root)
            self.CONNECTED = True
            self.Log(self.site_root)
        except Exception as ex:
            print('\nConnection to {} failed!'.format(self.site_root))
            print(ex)
            self.CONNECTED = False
            self.Log('ERROR: %s' % self.site_root)
            raise RuntimeError

    def Get(self, query):
        """
        Perform an http GET using:

         - *query* e.g.
         - *reply_mode* [default=''] this is the reply mode

        For example "/semanticSBML/annotate/search.xml?q=ATP"

        """
        data1 = None
        if self.CONNECTED:
            try:
                print(query)
                HTMLhead = {'User-Agent': self.USER_AGENT}
                self.conn.request("GET", query, headers=HTMLhead)
                r1 = self.conn.getresponse()
                print(r1.status, r1.reason)
                data1 = r1.read()
                self.Log('GET %s' % (query))
            except Exception as ex:
                print('\nFailure to GET: {}{}\n'.format(self.site_root, query))
                print(ex)
                self.Log('ERROR: %s%s' % (self.site_root, query))
                raise RuntimeError
        return data1

    def Close(self):
        """
        Close the currently active connection
        """
        if self.CONNECTED:
            self.conn.close()
            self.conn = None
            self.CONNECTED = False
            self.Log('%s - connection closed' % self.site_root)
            self.site_root = None


class MIRIAMTools(object):
    """
    Tools dealing with MIRIAM annotations
    """

    def MiriamURN2IdentifiersURL(self, urn):
        urn = urn.replace('urn:miriam:', '').split(':', 1)
        urn = 'http://identifiers.org/%s/%s' % (urn[0].strip(), urn[1].strip())
        print(urn)
        return urn


class SemanticSBML(RESTClient, MIRIAMTools):
    """
    REST client for connecting to SemanticSBML services

    """

    data = None
    item_re = re.compile('<item>.+?</item>')

    def __init__(self):
        RESTClient.__init__(self)

    def quickLookup(self, txt):
        """
        Do a quick lookpup for txt using SemanticSBML (connectic if required) and return results. Returns
        a list of identifiers.org id's in descending priority (as return)

         - *txt* the string to lookup

        """

        if not self.CONNECTED:
            self.Connect("www.semanticsbml.org")

        txt = txt.strip().replace(' ', '+')
        self.data = self.Get("/semanticSBML/annotate/search.xml?q={}".format(txt))
        self.data = self.parseXMLtoText(self.data)
        return self.data

    def viewDataInWebrowser(self, maxres=10):
        """
        Attempt to view #maxres results returned by SemanticSBML in the default browser

         - *maxres* default maximum number of results to display.

        """
        cntr = 0
        for u_ in self.data:
            cntr += 1
            try:
                webbrowser.open_new_tab(u_)
            except:
                print('ERROR in url {}'.format(u_))
            if cntr >= maxres:
                print('\nMaximum results reached, set \"maxres\" to increase')
                break

    def parseXMLtoText(self, xml):
        """
        Parse the xml output by quickLookup() into a list of URL

         - *xml* XML returns from SemanticSBML

        """
        return [
            i.replace('<item>', '').replace('</item>', '').strip()
            for i in re.findall(self.item_re, xml)
        ]

package ge.dchechelashvili.notesapp.main

import android.os.AsyncTask
import ge.dchechelashvili.notesapp.data.ItemsDatabase
import ge.dchechelashvili.notesapp.data.NotesDatabase
import ge.dchechelashvili.notesapp.data.entity.Item
import ge.dchechelashvili.notesapp.data.entity.Note

class MainInteractor(val presenter: IMainPresenter) {
    fun getNoteListFromDatabase(filterString: String) {
        GetNotesTask(presenter, filterString).execute()
    }

    fun getItemListFromDatabase(notes: List<Note>) {
        getItemsTask(presenter, notes).execute()
    }

    class GetNotesTask(val presenter: IMainPresenter, val filterString: String) : AsyncTask<Void, Void, List<Note>>() {
        override fun doInBackground(vararg params: Void?): List<Note> {
            val notesDao = NotesDatabase.getInstance().notesDao()
            val list = notesDao.getNotes(filterString)
            return list
        }

        override fun onPostExecute(result: List<Note>?) {
            super.onPostExecute(result)
            if (result != null) {
                presenter.onNoteListFetched(result)
            }
        }

    }

    class getItemsTask(val presenter: IMainPresenter, val notes: List<Note>) :
        AsyncTask<Void, Void, List<List<Item>>>() {
        override fun doInBackground(vararg params: Void?): List<List<Item>> {
            val itemsDao = ItemsDatabase.getInstance().itemsDao()
            val items = mutableListOf<List<Item>>()
            for (note in notes) {
                val list = itemsDao.getNoteItems(note.id)
                items.add(list)
            }
            return items
        }

        override fun onPostExecute(result: List<List<Item>>?) {
            super.onPostExecute(result)
            if (result != null) {
                presenter.onItemListFetched(result)
            }
        }

    }
}
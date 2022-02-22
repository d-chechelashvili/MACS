package ge.dchechelashvili.notesapp.details

import android.os.AsyncTask
import android.util.Log
import ge.dchechelashvili.notesapp.data.ItemsDatabase
import ge.dchechelashvili.notesapp.data.NotesDatabase
import ge.dchechelashvili.notesapp.data.entity.Item
import ge.dchechelashvili.notesapp.data.entity.Note
import ge.dchechelashvili.notesapp.main.IMainPresenter

class DetailsInteractor(val presenter: IDetailsPresenter) {

    fun saveInfoToDatabase(noteId: Int, note: Note, items: List<Item>) {
        SaveInfoTask(noteId, note, items, presenter).execute()
    }

    fun getNoteFromDatabase(id: Int) {
        GetNoteTask(presenter, id).execute()
    }

    fun getItemsFromDatabase(id: Int) {
        getNoteItemsTask(presenter, id).execute()
    }

    class getNoteItemsTask(val presenter: IDetailsPresenter, val id: Int) :
        AsyncTask<Void, Void, List<Item>>() {
        override fun doInBackground(vararg params: Void?): List<Item> {
            val itemsDao = ItemsDatabase.getInstance().itemsDao()
            val items = itemsDao.getNoteItems(id)
            return items
        }

        override fun onPostExecute(result: List<Item>?) {
            super.onPostExecute(result)
            if (result != null) {
                presenter.onNoteItemsFetched(result)
            }
        }

    }

    class GetNoteTask(val presenter: IDetailsPresenter, val id: Int) : AsyncTask<Void, Void, Note>() {

        override fun doInBackground(vararg params: Void?): Note{
            val notesDao = NotesDatabase.getInstance().notesDao()
            val note = notesDao.getNote(id)
            return note
        }

        override fun onPostExecute(result: Note?) {
            super.onPostExecute(result)
            if (result != null) {
                presenter.onNoteFetched(result)
            }
        }

    }


    class SaveInfoTask(val noteId: Int, val note: Note, val items: List<Item>, val presenter: IDetailsPresenter): AsyncTask<Void, Void, Void>(){
        override fun doInBackground(vararg params: Void?): Void? {
            val notesDao = NotesDatabase.getInstance().notesDao()
            val itemsDao = ItemsDatabase.getInstance().itemsDao()
            if (noteId != -1){
                Log.d("IDI", noteId.toString())
                notesDao.removeNote(noteId)
            }
            val noteId = notesDao.updateNote(note).toInt()
            itemsDao.removeItems(noteId)
            for (item in items){
                val newItem = Item(item.content, noteId, item.checked)
                itemsDao.updateItem(newItem)
            }
            return null
        }

        override fun onPostExecute(result: Void?) {
            super.onPostExecute(result)
            presenter.onInfoSaved()
        }

    }
}
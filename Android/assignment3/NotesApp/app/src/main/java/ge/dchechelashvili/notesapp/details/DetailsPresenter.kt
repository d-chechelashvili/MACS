package ge.dchechelashvili.notesapp.details

import android.util.Log
import ge.dchechelashvili.notesapp.data.entity.Item
import ge.dchechelashvili.notesapp.data.entity.Note
import java.util.*

class DetailsPresenter(val view: IDetailsView?): IDetailsPresenter {
    private val interactor = DetailsInteractor(this)

    override fun saveInfo(
        noteId: Int,
        name: String,
        pinned: Boolean,
        checkedItems: MutableList<Item>,
        uncheckedItems: MutableList<Item>
    ) {
        val note = Note(name, Calendar.getInstance().timeInMillis, pinned)
        val items = checkedItems + uncheckedItems
        interactor.saveInfoToDatabase(noteId, note, items)
    }

    override fun onInfoSaved() {
        view?.exit()
    }

    override fun getNote(id: Int) {
        interactor.getNoteFromDatabase(id)
    }

    override fun onNoteFetched(note: Note) {
        view?.showNote(note)
        interactor.getItemsFromDatabase(note.id)
    }

    override fun onNoteItemsFetched(items: List<Item>) {
        val checkedItems = mutableListOf<Item>()
        val uncheckedItems = mutableListOf<Item>()

        for (item in items){
            if(item.checked){
                checkedItems.add(item)
            }else{
                uncheckedItems.add(item)
            }
        }
        view?.showCheckedItems(checkedItems)
        view?.showUncheckedItems(uncheckedItems)
    }

}
package ge.dchechelashvili.notesapp.main

import ge.dchechelashvili.notesapp.data.entity.Item
import ge.dchechelashvili.notesapp.data.entity.Note

class MainPresenter(var view: IMainView?) : IMainPresenter {
    private val interactor = MainInteractor(this)
    private lateinit var notes: List<Note>

    fun getNotes(filterString: String) {
        interactor.getNoteListFromDatabase(filterString)
    }

    override fun onNoteListFetched(notes: List<Note>) {
        this.notes = notes
        interactor.getItemListFromDatabase(notes)
    }

    override fun onItemListFetched(items: List<List<Item>>) {
        val pinnedNotes = mutableListOf<Note>()
        val pinnedItems = mutableListOf<List<Item>>()
        val unpinnedNotes = mutableListOf<Note>()
        val unpinnedItems = mutableListOf<List<Item>>()

        for ((index, note) in notes.withIndex()) {
            val item = items[index]
            if (note.pinned) {
                pinnedNotes.add(note)
                pinnedItems.add(item)
            } else {
                unpinnedNotes.add(note)
                unpinnedItems.add(item)
            }
        }

        view?.showPinnedNotes(pinnedNotes, pinnedItems)
        view?.showUnpinnedNotes(unpinnedNotes, unpinnedItems)
    }


    fun detachView() {
        view = null
    }

}